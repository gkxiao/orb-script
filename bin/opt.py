#!/usr/bin/env python3
"""
Optimize molecular geometry with optional dihedral angle constraint using ORB v3 force field.

Version: 0.14
Purpose:
    - Read a molecule from XYZ file.
    - Optionally apply dihedral constraint (if --dihedral_indices provided).
    - If --dihedral_angle is not provided but indices are, use current dihedral angle as target.
    - Optimize geometry using ASE + ORB V3 force field.
    - Output optimized structure with relevant info in comment line:
        - With constraint: dihedral, actual_angle, energy, charge, spin, constraint_source
        - Without constraint: energy, charge, spin

Usage:
    # Unconstrained optimization:
    python opt_v0.14.py --input_file input.xyz --output_file output.xyz --charge 0.0 --spin 1.0

    # Freeze current dihedral angle:
    python opt_v0.14.py --input_file input.xyz --output_file output.xyz --charge 0.0 --spin 1.0 --dihedral_indices 2 3 4 12

    # Constrain to specific angle:
    python opt_v0.14.py --input_file input.xyz --output_file output.xyz --charge 0.0 --spin 1.0 --dihedral_indices 2 3 4 12 --dihedral_angle -90.0

Author: Gaokeng Xiao, Guangzhou Molcalx Ltd.
Date: 2025-04-05
"""

__version__ = "0.14"

import ase
from ase.build import bulk
from ase.io import read, write
from ase.optimize import BFGS
from orb_models.forcefield import pretrained
from orb_models.forcefield.calculator import ORBCalculator
from ase.io import Trajectory
from ase.constraints import FixInternals
from ase import Atoms
import numpy as np
import argparse

# ✅ 手动实现二面角计算，不依赖 ASE 内部函数
def calc_dihedral(p1, p2, p3, p4):
    """
    计算四个点 p1-p2-p3-p4 定义的二面角（单位：度）
    使用标准化学定义（右手法则）
    """
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)

    cos_phi = np.clip(np.dot(n1, n2), -1.0, 1.0)
    phi = np.arccos(cos_phi)

    sign = np.sign(np.dot(b2, np.cross(n1, n2)))
    if sign != 0:
        phi = phi * sign

    angle_deg = np.degrees(phi)
    return angle_deg

# --- 命令行参数解析 ---
parser = argparse.ArgumentParser(description="Optimize molecule with optional dihedral constraint.")
parser.add_argument("--input_file", type=str, required=True, help="Input XYZ file")
parser.add_argument("--output_file", type=str, required=True, help="Output XYZ file")
parser.add_argument("--charge", type=float, default=0.0, help="Total charge (default: 0.0)")
parser.add_argument("--spin", type=float, default=1.0, help="Spin multiplicity (2S+1, default: 1.0)")
parser.add_argument("--dihedral_indices", type=int, nargs=4, default=None,
                    help="Four 1-based atom indices defining dihedral (e.g., 2 3 4 12). If not provided, no constraint is applied.")
parser.add_argument("--dihedral_angle", type=float, nargs='?', default=None,
                    help="Target dihedral angle in degrees [0, 360]. Only used if --dihedral_indices is provided. If not given, current angle is used.")

args = parser.parse_args()

print(f"Running {__file__} version {__version__}")
print("-" * 60)

# 读入分子
atoms = read(args.input_file)
atoms.info["charge"] = args.charge
atoms.info["spin"] = args.spin

# 判断是否启用二面角约束
use_dihedral_constraint = args.dihedral_indices is not None

if use_dihedral_constraint:
    dihedral_indices_0based = [i - 1 for i in args.dihedral_indices]
    print(f"Constraining dihedral {args.dihedral_indices} (1-based) → {dihedral_indices_0based} (0-based)")

    # 自动决定目标角度
    if args.dihedral_angle is None:
        i, j, k, l = dihedral_indices_0based
        p = atoms.get_positions()
        current_angle = calc_dihedral(p[i], p[j], p[k], p[l])
        target_angle = current_angle
        print(f"--dihedral_angle not specified. Using current dihedral angle: {current_angle:.2f}°")
    else:
        target_angle = args.dihedral_angle
        print(f"Using user-specified target dihedral angle: {target_angle:.2f}°")

    # 创建并应用约束
    dihedral_constraint = [(target_angle, dihedral_indices_0based)]
    constraint = FixInternals(dihedrals_deg=dihedral_constraint)
    atoms.set_constraint(constraint)
else:
    print("No dihedral constraint applied. Performing unconstrained optimization.")
    atoms.set_constraint([])  # 清除任何已有约束（安全）

# 设置计算器
device = "cpu"
orbff = pretrained.orb_v3_conservative_inf_omat(
    device=device,
    precision="float32-high",
)
calc = ORBCalculator(orbff, device=device)
atoms.calc = calc

# 优化
print("Starting geometry optimization...")
dyn = BFGS(atoms)
dyn.run(fmax=0.05)

# 计算最终能量
final_energy = atoms.get_potential_energy()

# 仅在有约束时计算实际角度
actual_angle_deg = None
if use_dihedral_constraint:
    i, j, k, l = dihedral_indices_0based
    p = atoms.get_positions()
    actual_angle_deg = calc_dihedral(p[i], p[j], p[k], p[l])
    actual_angle_deg = actual_angle_deg % 360
    if actual_angle_deg < 0:
        actual_angle_deg += 360
    print(f"Actual dihedral angle: {actual_angle_deg:.2f}°")

print(f"Optimization finished.")
print(f"Final energy: {final_energy:.6f} eV")

# 创建干净原子对象
clean_atoms = Atoms(
    symbols=atoms.get_chemical_symbols(),
    positions=atoms.get_positions().copy(),
    cell=atoms.get_cell().copy(),
    pbc=atoms.get_pbc().copy()
)

# ✅ 动态构造 comment 行
if use_dihedral_constraint:
    dihedral_str = "-".join(map(str, args.dihedral_indices))
    angle_source = "fixed_current" if args.dihedral_angle is None else "user_specified"
    comment = f"dihedral={dihedral_str} actual_angle={actual_angle_deg:.2f} energy={final_energy:.6f} eV charge={args.charge} spin={args.spin} constraint_source={angle_source}"
else:
    comment = f"energy={final_energy:.6f} eV charge={args.charge} spin={args.spin}"

# 写入 XYZ 文件
with open(args.output_file, 'w') as f:
    f.write(f"{len(clean_atoms)}\n")
    f.write(f"{comment}\n")
    for atom in clean_atoms:
        f.write(f"{atom.symbol:2s} {atom.position[0]:15.8f} {atom.position[1]:15.8f} {atom.position[2]:15.8f}\n")

print(f"Optimized structure saved to {args.output_file} with comment: {comment}")
