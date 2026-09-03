#!/usr/bin/env bash
# =============================================================================
#  run_md_pipeline.sh
#
#  Runs the complete MD workflow for one solvated system, in order:
#
#     1. minimization        steepest descent, restrained
#     2. npt1 / nvt1  ...  npt7 / nvt7    7 equilibration cycles,
#                                         heating 100 K -> 303.15 K while the
#                                         position restraints are released
#     3. pre_production      50 ns unrestrained NPT
#     4. production          N independent replicas of 500 ns
#
#  Each step starts from the final structure (.gro) of the previous one; the
#  chain of steps is built automatically by the loops below.
#
# -----------------------------------------------------------------------------
#  HOW TO USE IT
#
#  Put this file next to your system files and run it from that directory:
#
#      cd Cripto-1_10D1/
#      ../scripts/run_md_pipeline.sh
#
#  To keep it running after you close the terminal:
#
#      nohup ../scripts/run_md_pipeline.sh > pipeline.log 2>&1 &
#
#  Everything is configurable from the command line, no need to edit the file:
#
#      NT=24 ../scripts/run_md_pipeline.sh                       # 24 threads
#      REPLICAS=5 ../scripts/run_md_pipeline.sh                  # 5 replicas
#      MDRUN_OPTS="-ntmpi 1 -ntomp 12 -nb gpu" ../scripts/run_md_pipeline.sh
#
#  If it crashes or you kill it, just run it again: steps that already finished
#  are skipped, and an interrupted step restarts from its GROMACS checkpoint.
# =============================================================================


# -----------------------------------------------------------------------------
#  Safety settings (three flags that make bash behave sensibly)
# -----------------------------------------------------------------------------
#   -e           stop the whole script as soon as any command fails, instead of
#                blindly continuing (if npt2 crashes, the script stops there
#                rather than launching npt3 on a file that was never written)
#   -u           error out if a variable is used but never defined (typos)
#   -o pipefail  a command pipeline fails if ANY part of it fails, not just the
#                last one
set -euo pipefail


# -----------------------------------------------------------------------------
#  CONFIGURATION
#
#  The "${VAR:-default}" syntax means: use VAR if it was given on the command
#  line, otherwise use the default written here.
# -----------------------------------------------------------------------------

GMX=${GMX:-gmx}                        # how GROMACS is called: gmx, gmx_mpi, ...
MDP_DIR=${MDP_DIR:-../MDP}             # folder containing all the .mdp files
TOP=${TOP:-topol.top}                  # topology of this system
NDX=${NDX:-index.ndx}                  # index file: MUST define SOLU and SOLV
START=${START:-step3_input.gro}        # solvated + ionized system from CHARMM-GUI.
                                       # Use the .gro, not the .pdb: the .gro
                                       # carries the periodic box vectors on its
                                       # last line, which grompp needs. Serves as
                                       # both coordinates (-c) and restraint
                                       # reference (-r) for the minimization.

CYCLES=${CYCLES:-7}                    # number of NPT/NVT equilibration cycles
REPLICAS=${REPLICAS:-3}                # number of independent production runs
PROD_NSTEPS=${PROD_NSTEPS:-250000000}  # 250e6 steps x 2 fs = 500 ns per replica

NT=${NT:-48}                           # CPU threads for mdrun
MDRUN_OPTS=${MDRUN_OPTS:--nt ${NT}}    # any extra mdrun flags (GPU, pinning...)

LOG_DIR=${LOG_DIR:-logs}               # all grompp/mdrun output goes here
MAXWARN=${MAXWARN:-0}                  # grompp warnings tolerated (keep at 0!)


# -----------------------------------------------------------------------------
#  Small helper functions
# -----------------------------------------------------------------------------

# Print a timestamped message, so the log tells you when each step started.
say() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# Print an error and stop the script.
die() {
    echo "" >&2
    echo "ERROR: $*" >&2
    exit 1
}

# Check that a file exists before we try to use it.
check_file() {
    [[ -f "$1" ]] || die "missing file: $1  (are you in the right directory?)"
}

# Counters used only to print "step 5 of 18" while running.
STEP_NUMBER=0
TOTAL_STEPS=$(( 1 + 2 * CYCLES + 1 + REPLICAS ))


# -----------------------------------------------------------------------------
#  run_step - the core of the script
#
#  It runs one simulation step: grompp (build the run input) followed by
#  mdrun (run it). It is called once per step of the workflow.
#
#      run_step <name> <input structure> [mdp file]
#
#  <name> is the base name used for every file that step produces:
#  name.tpr, name.gro, name.xtc, name.log, name.cpt ...
#  If no .mdp is given, it looks for MDP_DIR/<name>.mdp.
# -----------------------------------------------------------------------------
run_step() {
    local name="$1"                                   # e.g. npt3
    local input="$2"                                  # e.g. nvt2.gro
    local mdp="${3:-${MDP_DIR}/$(basename "${name}").mdp}"

    # "name" can contain a subfolder (production/rep1/production); for log file
    # names we replace the slashes with underscores.
    local tag="${name//\//_}"

    STEP_NUMBER=$(( STEP_NUMBER + 1 ))

    # --- Already done? Then there is nothing to do. -------------------------
    # This is what makes the script restartable: the .gro file is only written
    # when a step has completed successfully.
    if [[ -f "${name}.gro" ]]; then
        say "[${STEP_NUMBER}/${TOTAL_STEPS}] ${name}: already finished, skipping"
        return 0
    fi

    mkdir -p "$(dirname "${name}")" "${LOG_DIR}"
    check_file "${mdp}"
    check_file "${input}"

    say "[${STEP_NUMBER}/${TOTAL_STEPS}] ${name}: preparing (grompp)"

    # --- grompp: combine structure + topology + parameters into a .tpr ------
    # -c : starting coordinates
    # -r : reference coordinates for the position restraints (same file here)
    # We skip this if the .tpr is already there from a previous attempt.
    if [[ ! -f "${name}.tpr" ]]; then
        ${GMX} grompp -f "${mdp}" \
                      -o "${name}.tpr" \
                      -c "${input}" \
                      -r "${input}" \
                      -p "${TOP}" \
                      -n "${NDX}" \
                      -maxwarn "${MAXWARN}" \
                      > "${LOG_DIR}/grompp_${tag}.log" 2>&1 \
            || {
                echo "--- last lines of ${LOG_DIR}/grompp_${tag}.log ---" >&2
                tail -n 30 "${LOG_DIR}/grompp_${tag}.log" >&2
                die "grompp failed at step '${name}'"
            }
    fi

    # --- mdrun: run the simulation ------------------------------------------
    # If a checkpoint (.cpt) exists, the previous attempt was interrupted:
    # continue from where it stopped instead of starting over.
    local restart=""
    if [[ -f "${name}.cpt" ]]; then
        restart="-cpi ${name}.cpt"
        say "           found a checkpoint, resuming from it"
    fi

    say "[${STEP_NUMBER}/${TOTAL_STEPS}] ${name}: running (mdrun)"

    # ${restart} and ${MDRUN_OPTS} are intentionally NOT quoted: we want bash
    # to split them into separate arguments (e.g. "-nt" and "48").
    ${GMX} mdrun -deffnm "${name}" ${restart} ${MDRUN_OPTS} \
        > "${LOG_DIR}/mdrun_${tag}.log" 2>&1 \
        || {
            echo "--- last lines of ${LOG_DIR}/mdrun_${tag}.log ---" >&2
            tail -n 30 "${LOG_DIR}/mdrun_${tag}.log" >&2
            die "mdrun failed at step '${name}'"
        }

    say "[${STEP_NUMBER}/${TOTAL_STEPS}] ${name}: done"
}


# =============================================================================
#  START HERE
# =============================================================================

# --- Check everything is in place before starting a multi-day job ------------
command -v "${GMX}" >/dev/null || die "'${GMX}' not found. Did you load the GROMACS module?"
[[ -d "${MDP_DIR}" ]]          || die "mdp folder not found: ${MDP_DIR}"
check_file "${TOP}"
check_file "${NDX}"
check_file "${START}"

mkdir -p "${LOG_DIR}"

echo "============================================================"
say  "starting MD pipeline"
echo "  working directory : $(pwd)"
echo "  GROMACS           : $(${GMX} --version | awk '/^GROMACS version/ {print $3}')"
echo "  mdp files         : ${MDP_DIR}"
echo "  starting structure: ${START}"
echo "  mdrun options     : ${MDRUN_OPTS}"
echo "  equilib. cycles   : ${CYCLES}"
echo "  production        : ${REPLICAS} replicas x ${PROD_NSTEPS} steps"
echo "  logs              : ${LOG_DIR}/"
echo "============================================================"


# -----------------------------------------------------------------------------
#  1. Energy minimization
#
#  A single steepest-descent run (minimization.mdp), starting from the
#  CHARMM-GUI system. The same file is used as coordinates (-c) and as the
#  reference for the position restraints (-r), since nothing has moved yet.
#  It produces minimization.gro, the input of the first equilibration cycle.
# -----------------------------------------------------------------------------
run_step minimization "${START}"


# -----------------------------------------------------------------------------
#  2. Equilibration: CYCLES x (NPT then NVT)
#
#  "previous" always holds the structure the next step has to start from:
#  minimization -> npt1 -> nvt1 -> npt2 -> nvt2 -> ... -> nvt7
# -----------------------------------------------------------------------------
previous="minimization.gro"

for i in $(seq 1 "${CYCLES}"); do
    echo ""
    say "--- equilibration cycle ${i} of ${CYCLES} ---"

    run_step "npt${i}" "${previous}"      # constant pressure
    run_step "nvt${i}" "npt${i}.gro"      # constant volume

    previous="nvt${i}.gro"
done


# -----------------------------------------------------------------------------
#  3. Pre-production: 50 ns, no restraints
# -----------------------------------------------------------------------------
echo ""
say "--- pre-production ---"
run_step pre_production "${previous}"


# -----------------------------------------------------------------------------
#  4. Production replicas
#
#  production_replicas.mdp is distributed with "nsteps = -1" (the original runs
#  were extended manually). Here we write a copy with the real length, so each
#  replica stops by itself after PROD_NSTEPS steps.
# -----------------------------------------------------------------------------
echo ""
say "--- production ---"

production_mdp="${MDP_DIR}/production_replicas.mdp"

if [[ -n "${PROD_NSTEPS}" ]]; then
    mkdir -p production
    production_mdp="production/production.mdp"
    sed -E "s/^[[:space:]]*nsteps[[:space:]]*=.*/nsteps = ${PROD_NSTEPS}/" \
        "${MDP_DIR}/production_replicas.mdp" > "${production_mdp}"
    say "production length set to ${PROD_NSTEPS} steps in ${production_mdp}"
fi

# Every replica starts from the same pre-production structure but gets its own
# grompp call. Since the mdp has "gen-vel = yes" and "gen-seed = -1", GROMACS
# draws a new random velocity distribution each time - that is exactly what
# makes the three replicas independent.
for r in $(seq 1 "${REPLICAS}"); do
    run_step "production/rep${r}/production" "pre_production.gro" "${production_mdp}"
done


echo ""
echo "============================================================"
say  "pipeline finished successfully"
echo "  trajectories: production/rep*/production.xtc"
echo "============================================================"
