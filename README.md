## Introduction
Orb-v3 <sup>1,2</sup> is a next-generation universal machine learning interatomic potential (MLIP) developed by Orbital Materials, achieving significant breakthroughs on the Pareto frontier of performance-speed-memory trade-offs. This model series includes multiple variants covering different configurations from high-precision conservative types to high-efficiency sparse graph non-conservative types, making it suitable for a wide range of computational chemistry tasks.

The scripts are built upon the Orb-V3 force field model.

## Tutorial
<img src="http://blog.molcalx.com.cn/wp-content/uploads/2025/09/twisted_dihedral_angle-2-3-4-12-768x494.png">
<p style="text-align:center;">Figure 1. Example molecule. Highlighted dihedral: 2-3-4-12 (1-based indices)</p>
➤ regular optimization：
```bash
python opt.py --input_file lig_twisted.xyz --output_file opted.xyz --charge 0 --spin 1
```

