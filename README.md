## Introduction
Orb-v3 <sup>1,2</sup> is a next-generation universal machine learning interatomic potential (MLIP) developed by Orbital Materials, achieving significant breakthroughs on the Pareto frontier of performance-speed-memory trade-offs. This model series includes multiple variants covering different configurations from high-precision conservative types to high-efficiency sparse graph non-conservative types, making it suitable for a wide range of computational chemistry tasks.

The scripts are built upon the Orb-V3 force field model.

## Example Molecule
<img src="http://blog.molcalx.com.cn/wp-content/uploads/2025/09/twisted_dihedral_angle-2-3-4-12-768x494.png" style="width:350px; height:225px;">
<p style="text-align:center;">Figure 1. Example molecule. Highlighted dihedral: 2-3-4-12 (1-based indices)</p>

<p>The test molecule (Figure 1) is the ligand extracted from the BTK co-crystal structure PDB 4ZLZ. Its dihedral angle 2-3-4-12 is rotated from the bioactive conformation of -73.8° to -139.1°.</p>
<p>The files for these angles are provided in the data directory.</p>
<ul>
<li>lig_twisted.xyz:  dihedral angle 2-3-4-12 = -139.1°</li>
<li>lig_xtal.xyz:  dihedral angle 2-3-4-12 = -73.8°</li>
</ul>

## Optimization tutorial
➤ regular optimization：
```bash
python opt.py --input_file lig_twisted.xyz --output_file opt.xyz --charge 0 --spin 1
```
➤ Optimization with specific dihedral angles frozen
```bash
python opt.py --input_file opt.xyz --output_file opt_frozen.xyz --charge 0 --spin 1 --dihedral_indices 2 3 4 12
```
➤ Optimization with specific dihedral angle constraints, such as -90 degrees
```bash
python opt.py --input_file lig_twisted.xyz --output_file opt_-90.xyz --charge 0 --spin 1 --dihedral_indices 2 3 4 12 --dihedral_angle -90
```
## Cluster tutorial
➤ Geometric optimization of the conformational ensemble was performed using Orb V3, followed by clustering with [CREST (V3.01)](https://github.com/crest-lab/crest)
using energy windows = 10 kcal/mol and RMSD threshold = 0.125 Å to remove duplicate conformations.
```bash
cluster.py confs.sdf --charge 0
```
<p>where the parameter charge refers to the net formal charge of the molecule.</p>
## Reference
<ol>
    <li>Rhodes, B. et al. (2025) “Orb-v3: atomistic simulation at scale.” Available at: http://arxiv.org/abs/2504.06231.</li>
    <li>https://github.com/orbital-materials/orb-models</li>
</ol>
