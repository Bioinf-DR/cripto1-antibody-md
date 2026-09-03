# MD input files — h10D1/Cripto-1 and h1B4/Cripto-1 complexes

Input files needed to reproduce the all-atom molecular dynamics simulations reported in:

> *Distinct molecular recognition of Cripto-1 by two humanized antibodies: a combined experimental
> and computational investigation.* Oliver A., Migliorini A., Colletta M., Iaccarino E., Barisciano G.,
> Sanna R., Ronca R., Gracia Carmona O., Ruvo M., Raimondo D., Sandomenico A.
> *International Journal of Biological Macromolecules* (2026). DOI: `<TO FILL>`

This deposit contains the **starting coordinates, topologies and run parameters** for the two
simulated systems, plus the driver scripts used for minimization, equilibration and production.
Trajectories are not included; they are available from the corresponding author on reasonable
request.

---

## 1. Systems

| Directory | System | Contents |
|---|---|---|
| `h10D1_Cripto1/` | Fab h10D1 in complex with Cripto-1 (residues 31–161) | starting structure, topology, force-field files, index file |
| `h1B4_Cripto1/`  | Fab h1B4 in complex with Cripto-1 (residues 31–161)  | starting structure, topology, force-field files, index file |

Both starting structures are the HADDOCK 2.4 top-scoring docking poses, solvated and ionized with
the CHARMM-GUI *Solution Builder*.

**Preparation summary**

- Antigen: full-length Cripto-1 (UniProt P13385) predicted with AlphaFold 2.3; the first 30 residues
  (signal peptide) were removed, leaving the mature soluble form 31–161. The N-terminal residue is
  capped with an acetyl group.
- Antibodies: Fv/Fab models of h10D1 and h1B4 built with RosettaAntibody (ROSIE) and refined at
  CDR-H3 with RFdiffusion.
- Force field: **CHARMM36m** for the protein, **TIP3P** water.
- Box: rectangular, periodic, 10 Å solute–box buffer.
- Ions: 0.15 M NaCl, system neutral.
- Protonation states from PropKa at pH 7. Histidines protonated on Nδ: H33, H74, H108, H135, H157
  (Cripto-1); H38 (h10D1 heavy chain); H107 (h10D1 light chain). Cripto-1 **H120 is ε-protonated**.
- Disulfides: Cripto-1 C115–C133, C128–C149, C131–C140; both Fabs C23–C104 in heavy and light chain.

Software: **GROMACS 2022.3**.

---

## 2. Directory layout

```
.
├── README.md
├── mdp/
│   ├── minimization.mdp
│   ├── npt1.mdp … npt7.mdp
│   ├── nvt1.mdp … nvt7.mdp
│   ├── pre_production.mdp
│   └── production_replicas.mdp
├── scripts/
│   ├── execute_minimization.sh
│   └── execute_equilibration_preproduction.sh
├── h10D1_Cripto1/
│   ├── step3_input.pdb     # solvated, ionized system (CHARMM-GUI output)
│   ├── topol.top
│   ├── toppar/             # CHARMM36m itp files written by CHARMM-GUI
│   ├── *.itp               # chain topologies and position-restraint files
│   └── index.ndx           # must define the SOLU and SOLV groups
└── h1B4_Cripto1/
    └── (same layout)
```
---

## 3. Simulation protocol

### 3.1 Energy minimization — `minimization.mdp`

Steepest descent, 5000 steps, `emtol = 1000 kJ mol⁻¹ nm⁻¹`, with harmonic position restraints of
2000 kJ mol⁻¹ nm⁻² on backbone atoms and 200 kJ mol⁻¹ nm⁻² on side chains of both Cripto-1 and the
antibody.

### 3.2 Equilibration — seven NPT/NVT cycles

Temperature is ramped from 100 K to 303.15 K while position restraints are progressively released.
Temperature coupling: velocity-rescale thermostat, `tau_t = 1.0 ps`. Pressure coupling (NPT only):
stochastic C-rescale barostat, isotropic, `tau_p = 2.0 ps`, `ref_p = 1.0 bar`.

| Cycle | `.mdp` files | *T* (K) | Backbone restraint (kJ mol⁻¹ nm⁻²) | Side-chain restraint | NPT length | NVT length |
|---|---|---|---|---|---|---|
| 1 | `npt1` / `nvt1` | 100 | 2000 | 200 | 1 ns | 1 ns |
| 2 | `npt2` / `nvt2` | 150 | 1000 | 100 | 1 ns | 1 ns |
| 3 | `npt3` / `nvt3` | 200 | 500 | 50 | 1 ns | 1 ns |
| 4 | `npt4` / `nvt4` | 250 | 250 | 25 | 1 ns | 1 ns |
| 5 | `npt5` / `nvt5` | 280 | 150 | 10 | 1 ns | 1 ns |
| 6 | `npt6` / `nvt6` | 303.15 | 50 | 0 | 1 ns | 1 ns |
| 7 | `npt7` / `nvt7` | 303.15 | none | none | 2 ns | 2 ns |

Velocities are generated once at the start (`gen-vel = yes`, `gen-temp = 100` in `npt1.mdp`);
all later steps continue from the previous one (`continuation = yes`).

Common settings across all equilibration and production steps: `dt = 2 fs`; `cutoff-scheme = Verlet`
with `nstlist = 20`; van der Waals cut-off 1.2 nm with force-switch from 1.0 nm; PME electrostatics
with a 1.2 nm real-space cut-off; hydrogen bonds constrained with LINCS (`lincs_order = 4`,
`lincs_iter = 1`); 3D periodic boundary conditions.

### 3.3 Pre-production — `pre_production.mdp`

50 ns unrestrained NPT at 303.15 K, coordinates saved every 10 ps.

### 3.4 Production — `production_replicas.mdp`

Unrestrained NPT at 303.15 K, coordinates saved every 10 ps. **Three independent replicas of 500 ns
were run per system** (1.5 µs per system, 3 µs total). Replicas are made independent by
`gen-vel = yes` with `gen-seed = -1`, i.e. each replica draws a fresh set of velocities from a
Maxwell distribution.

`nsteps` is set to `-1` in the deposited file, since the runs were extended with `mdrun -maxh` on the
local queue. To reproduce a 500 ns replica, set:

```
nsteps = 250000000    ; 250000000 * 2 fs = 500 ns
```

---

## 4. How to rerun

From inside one of the system directories, with the `mdp/` files in place:

```bash
# 1. Minimization
gmx grompp -f minimization.mdp -o minimization.tpr \
           -c step3_input.pdb -r step3_input.pdb -p topol.top -n index.ndx
gmx mdrun -v -deffnm minimization

# 2. Equilibration + pre-production (7 NPT/NVT cycles, then 50 ns)
bash execute_equilibration_preproduction.sh

# 3. Production — repeat three times in separate directories
gmx grompp -f production_replicas.mdp -o production_replicas.tpr \
           -c pre_production.gro -r pre_production.gro -p topol.top -n index.ndx
gmx mdrun -v -deffnm production_replicas
```

`execute_equilibration_preproduction.sh` chains the steps in the order
`npt1 → nvt1 → npt2 → nvt2 → … → npt7 → nvt7 → pre_production`, each `grompp` taking the `.gro` of
the previous step as both `-c` and `-r`. Adjust `-nt 48` in the script to the core count of your
machine (or replace it with your scheduler's submission syntax).

Approximate wall time: `<TO FILL — e.g. ~X ns/day on N cores / GPU model>`.

---

## 5. Notes and known inconsistencies

Please read before rerunning:

1. **Do not deposit `execute_minimization.sh`.** That script loops over
   `minimization_2.mdp … minimization_21.mdp`, i.e. a 20-step staged minimization that was not the
   protocol used here: minimization was a single steepest-descent run with `minimization.mdp`,
   starting from `step3_input.pdb`.
2. **`npt4.mdp` sets `POSRES_FC_BB = 2500.0`** while the surrounding steps follow the sequence
   2000 → 1000 → 500 → **250** → 150 → 50, and `nvt4.mdp` uses 250.0. This looks like a typo that
   briefly *increased* backbone restraints at cycle 4. It does not invalidate the equilibration —
   the restraint is fully released by cycle 7 — but the manuscript states restraints were
   "progressively reduced across the cycles", so either fix the value or note the deviation.
3. **`nvt2.mdp` runs 0.5 ns**, not 1 ns like the other cycles. Total equilibration is therefore
   15.5 ns + 50 ns pre-production = **65.5 ns**, not the 66 ns quoted in the Methods section.
4. `pre_production.mdp` has `gen-vel = yes` without `gen-temp`, so velocities are regenerated at the
   GROMACS default (300 K) rather than 303.15 K. Harmless over a 50 ns run, but worth being aware of
   if you are comparing energies from the first few ps.
5. The position-restraint macros `-DPOSRES_FC_BB` and `-DPOSRES_FC_SC` require the CHARMM-GUI-style
   `posre*.itp` files; these must be present in the topology directory or `grompp` will silently
   apply no restraints.

---

## 6. Contact

`<TO FILL — corresponding author name and email>`

Licence: `<TO FILL — e.g. CC BY 4.0>`
