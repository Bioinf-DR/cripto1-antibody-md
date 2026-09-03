# MD input files — Cripto-1 in complex with the humanized antibodies 10D1 and 1B4

Input files needed to reproduce the all-atom molecular dynamics simulations reported in:

> *Distinct molecular recognition of Cripto-1 by two humanized antibodies: a combined experimental
> and computational investigation.* Oliver A., Migliorini A., Colletta M., Iaccarino E., Barisciano G.,
> Sanna R., Ronca R., Gracia Carmona O., Ruvo M., Raimondo D., Sandomenico A.
> *International Journal of Biological Macromolecules* (2026). DOI: `<TO FILL>`

This repository contains the **starting coordinates, topologies, run parameters and driver script**
for the two simulated systems. Trajectories are not included; they are available from the
corresponding author on reasonable request.

---

## 1. Systems

| Directory | System |
|---|---|
| `Cripto-1_10D1/` | h10D1 in complex with Cripto-1 (residues 31–161) |
| `Cripto-1_1B4/`  | h1B4 in complex with Cripto-1 (residues 31–161)  |

Both starting structures are the HADDOCK 2.4 top-scoring docking poses, solvated and ionized with
the CHARMM-GUI *Solution Builder*.

**Preparation summary**

- Antigen: full-length Cripto-1 (UniProt P13385) predicted with AlphaFold 2.3; the first 30 residues
  (signal peptide) were removed, leaving the mature soluble form 31–161. The N-terminal residue is
  capped with an acetyl group.
- Antibodies: models of h10D1 and h1B4 built with RosettaAntibody (ROSIE) and refined at CDR-H3 with
  RFdiffusion. Each simulated system contains three protein chains — Cripto-1 plus the antibody
  heavy and light chains — for a total of 360 residues.
  `<TO CONFIRM: 360 residues total corresponds to an Fv construct (VH + VL), not a Fab. Align this
  section and Methods 2.10 with whichever is correct.>`
- Force field: **CHARMM36m** for the protein, **TIP3P** water.
- Box: cubic, periodic, 11.6 nm per side (10 Å solute–box buffer).
- Ions: 0.15 M NaCl, system neutral.
- Protonation states from PropKa at pH 7. Histidines protonated on Nδ: H33, H74, H108, H135, H157
  (Cripto-1); H38 (h10D1 heavy chain); H107 (h10D1 light chain). Cripto-1 **H120 is ε-protonated**.
- Disulfides: Cripto-1 C115–C133, C128–C149, C131–C140; both antibodies C23–C104 in heavy and light
  chain.

Software used for the published runs: **GROMACS 2022.3**.

---

## 2. Repository layout

```
.
├── README.md
├── MDP/                          # run parameters, shared by both systems
│   ├── minimization.mdp
│   ├── npt1.mdp … npt7.mdp
│   ├── nvt1.mdp … nvt7.mdp
│   ├── pre_production.mdp
│   └── production_replicas.mdp
├── scripts/
│   └── run_md_pipeline.sh        # runs the whole workflow for one system
├── Cripto-1_10D1/
│   ├── step3_input.gro           # solvated, ionized system (CHARMM-GUI output)
│   ├── topol.top                 # includes everything in toppar/
│   ├── index.ndx                 # defines the SOLU and SOLV groups
│   └── toppar/                   # CHARMM36m parameters and chain topologies
│       ├── forcefield.itp
│       ├── PROA.itp              # `<TO CONFIRM: which chain is which>`
│       ├── PROB.itp
│       ├── PROC.itp
│       ├── TIP3.itp              # water
│       ├── SOD.itp               # Na+
│       └── CLA.itp               # Cl-
└── Cripto-1_1B4/
    └── (same layout)
```

> `topol.top` is only a wrapper: it `#include`s the files in `toppar/`, so the whole folder is
> required and must stay with its topology. The position-restraint blocks driven by
> `-DPOSRES_FC_BB` and `-DPOSRES_FC_SC` are defined inside the `PRO*.itp` chain topologies.
>
> Before a long run, check that everything resolves by building a test `.tpr` inside each system
> directory:
>
> ```bash
> cd Cripto-1_10D1
> gmx grompp -f ../MDP/minimization.mdp -o test.tpr \
>            -c step3_input.gro -r step3_input.gro -p topol.top -n index.ndx
> rm -f test.tpr mdout.mdp
> ```
>
> The system contains 146857 atoms: 5524 protein, 133 Na⁺, 137 Cl⁻ and 47021 waters. If `grompp`
> reports a different total, the structure and the topology do not belong to the same build.

---

## 3. Simulation protocol

### 3.1 Energy minimization — `minimization.mdp`

A single steepest-descent run, 5000 steps, `emtol = 1000 kJ mol⁻¹ nm⁻¹`, starting from
`step3_input.gro`, with harmonic position restraints of 2000 kJ mol⁻¹ nm⁻² on backbone atoms and
200 kJ mol⁻¹ nm⁻² on side chains of both Cripto-1 and the antibody.

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
`gen-vel = yes` with `gen-seed = -1`: each one draws a fresh set of velocities from a Maxwell
distribution.

`nsteps` is set to `-1` in the deposited file, since the runs were extended on the local queue.
A 500 ns replica corresponds to `nsteps = 250000000` (250 × 10⁶ × 2 fs); `run_md_pipeline.sh` sets
this automatically.

---

## 4. How to rerun

### With the pipeline script

From inside one of the system directories:

```bash
cd Cripto-1_10D1/
../scripts/run_md_pipeline.sh
```

This runs minimization → 7 NPT/NVT cycles → 50 ns pre-production → 3 production replicas, in order,
stopping if any step fails. To leave it running after closing the terminal:

```bash
nohup ../scripts/run_md_pipeline.sh > pipeline.log 2>&1 &
```

Settings are passed as environment variables, so the script does not need to be edited:

```bash
NT=24 ../scripts/run_md_pipeline.sh                                     # 24 CPU threads
REPLICAS=5 ../scripts/run_md_pipeline.sh                                # 5 replicas
GMX=gmx ../scripts/run_md_pipeline.sh                                 # double-precision build
MDRUN_OPTS="-ntmpi 1 -ntomp 12 -nb gpu" ../scripts/run_md_pipeline.sh   # GPU run
```

The script is restartable: finished steps are skipped and an interrupted step resumes from its
GROMACS checkpoint, so it is safe to launch again after a walltime kill.

Outputs are written in the system directory: `minimization.gro`, `npt*.gro`, `nvt*.gro`,
`pre_production.gro`, and the trajectories in `production/rep{1,2,3}/production.xtc`. All `grompp`
and `mdrun` output goes to `logs/`.

### Manually

Each step is a `grompp` followed by an `mdrun`, chained through the `.gro` files:

```bash
gmx grompp -f ../MDP/minimization.mdp -o minimization.tpr \
           -c step3_input.gro -r step3_input.gro -p topol.top -n index.ndx
gmx mdrun -deffnm minimization

gmx grompp -f ../MDP/npt1.mdp -o npt1.tpr \
           -c minimization.gro -r minimization.gro -p topol.top -n index.ndx
gmx mdrun -deffnm npt1
# … npt1 → nvt1 → npt2 → nvt2 → … → npt7 → nvt7 → pre_production
```
---

## 5. Notes

1. Total equilibration is 16 ns (six 1 ns + 1 ns cycles, then one 2 ns + 2 ns cycle), followed by
   50 ns of pre-production: 66 ns before the production runs begin.
2. `pre_production.mdp` has `gen-vel = yes` without `gen-temp`, so velocities are regenerated at the
   GROMACS default (300 K) rather than 303.15 K. This is inconsequential over a 50 ns run, but worth
   knowing if you compare energies from the first few picoseconds.
3. Start from `step3_input.gro` rather than a PDB: the `.gro` carries the periodic box vectors on its
   last line, which `grompp` requires.
4. The files were prepared with GROMACS 2022.3 and also build cleanly with 2025.2.

---

## 6. Contact and licence

`<TO FILL — corresponding author name and email>`

Licence: `<TO FILL — e.g. CC BY 4.0>`
