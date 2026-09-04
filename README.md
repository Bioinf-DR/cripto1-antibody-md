# MD input files — Cripto-1 in complex with the humanized antibodies 10D1 and 1B4

Input files needed to reproduce the all-atom molecular dynamics simulations reported in:

corresponding author on reasonable request.

---

## 1. Systems

| Directory | System |
|---|---|
| `Cripto-1_10D1/` | h10D1 in complex with Cripto-1 (residues 31–161) |
| `Cripto-1_1B4/`  | h1B4 in complex with Cripto-1 (residues 31–161)  |

Both starting structures are the HADDOCK 2.4 top-scoring docking poses, solvated and ionized with
the CHARMM-GUI *Solution Builder*.

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
│       ├── PROA.itp              
│       ├── PROB.itp
│       ├── PROC.itp
│       ├── TIP3.itp              # water
│       ├── SOD.itp               # Na+
│       └── CLA.itp               # Cl-
└── Cripto-1_1B4/
    └── (same layout)
```

---

## 3. How to rerun

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
3. Start from `step3_input.gro` rather than a PDB: the `.gro` carries the periodic box vectors on its
   last line, which `grompp` requires.

---

## 6. Contacts 

Arianna Migliorini: arianna.migliorini@uniroma1.it 
Miria Colletta: colletta.2004499@studenti.uniroma1.it

