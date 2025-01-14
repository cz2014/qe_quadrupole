# Quadrupole Patch for Quantum Espresso

A patch for [Quantum ESPRESSO](https://www.quantum-espresso.org/) (QE) that enables the computation of the quadrupole tensor from `ph.x`, which is not available in the original QE.

## Overview

This repository contains modified source files and example scripts to compute the quadrupole tensor for 3D or 2D materials. The main workflow involves:
1. Generating the **q-dependent charge densities** using the patched QE 6.7.
2. Processing the resulting data with `dyn2quad.py` to obtain a `quadrupole.fmt` file.
3. Incorporating `quadrupole.fmt` into subsequent DFPT and EPW calculations (with QE 6.8), and for any other use.

> **Note:** The current implementation assumes the material has 3-fold symmetry. If your material does not have 3-fold symmetry, disable the symmetry correction in `dyn2quad.py`.

---

## Usage

Below is a step-by-step guide on how to use the modified code and the accompanying scripts.

### 1. Generate q-Dependent Charge Densities (Using QE 6.7)

1. **Patch and compile QE 6.7** to include the quadrupole-extraction capability.  
2. Prepare an SCF input file (`scf.in`) with `nosym = .true.` to avoid computating with phonon pattern symmetry.  
3. Run `ph.x` at a **set of q points**. A full set includes 16 q points, but a minimum set of 7 q points \([q0, qx1, qy1, qz1, qxy1, qxz1, qyz1]\) is sufficient in most cases.  
4. Ensure your **k-point grid is sufficiently dense**, as the quadrupole tensor is sensitive to k-point sampling.

> **Note:** You can quickly check if your k-point grid is adequate by comparing total energies at different grids. For example, a 12×12×1 grid vs. an 18×18×1 grid for Sb can yield a 30% difference.

### 2. Collect Data and Generate `quadrupole.fmt`

1. Use the provided `dyn2quad.py` script to parse the phonon outputs (from step 1) and generate the file `quadrupole.fmt`.  
2. Verify that the **raw data points** and **fitted data** are consistent, by checking the generated figures from `dyn2quad.py` script. Large discrepancies might indicate issues such as forgetting to set `nosym = .true.`.

### 3. Use `quadrupole.fmt` in DFPT and EPW Calculations (Using QE 6.8)

1. Place `quadrupole.fmt` in the folders where `ph.x` and `epw.x` will be executed.  
2. In the DFPT calculation, set `assume_isolated = '2D'` and ensure the 2D material is centered around `z=0`.  
3. Use the **modified `q2r.x`** (provided in this repository) to obtain the force-constant file for EPW. The presence of `quadrupole.fmt` is required at this stage.  
4. For EPW runs, make sure `quadrupole.fmt` is also present so that quadrupole contributions are accurately accounted for.

---

## Examples

Two example workflows (`ex1_mos2/` and `ex2_inse/`) are included in this repository to demonstrate:
- How to structure input files
- How to run `scf.x` and `ph.x` at multiple q points
- How to use `dyn2quad.py` and integrate the resulting `quadrupole.fmt` into your DFPT and EPW calculations


---

## References

- Quantum ESPRESSO: [https://www.quantum-espresso.org/](https://www.quantum-espresso.org/)  
- Paper published: [https://doi.org/10.1103/PhysRevB.106.115423](https://doi.org/10.1103/PhysRevB.106.115423)
