#!/usr/bin/env bash
#SBATCH --job-name=Re100_T100
#SBATCH --partition=hpg-turin
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --gpus-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4-00:00:00
#SBATCH --output=gpu_%j.out
#SBATCH --error=gpu_%j.err

set -euo pipefail
# sbatch --chdir=<case> sets the actual working directory.
# Do NOT use $0 here: Slurm may execute a spooled copy of this script.
test -f np3d.par || { echo "Expected to run in its case directory (sbatch --chdir=...)" >&2; exit 11; }
exec 9>".stage5a_psweep_gpu.lock"
flock -n 9 || { echo "Case directory already locked" >&2; exit 1; }
test -s np3d.re2 && test -s physalis_cage_nodes.dat
module load cuda/13.2.1 gcc/14.2.0 openmpi/5.0.10 cmake/3.30.5
# Reuse previously built patched Stage5C dynamic core; never modify it.
export NEKRS_HOME="${NEKRS_PHYSALIS_CORE_HOME:?Set NEKRS_PHYSALIS_CORE_HOME to the validated patched NekRS install}"
if [[ ! -f "$NEKRS_HOME/.physalis_stage5c_h2_dynamic_core" || ! -x "$NEKRS_HOME/bin/nekrs" ]]; then
    echo "FATAL: requires previously built Stage5C dynamic core; export NEKRS_PHYSALIS_CORE_HOME to point at an existing solver-install" >&2
    exit 10
fi
export PATH="$NEKRS_HOME/bin:$PATH"
export NEKRS_PHYSALIS_CAGE_FILE="$PWD/physalis_cage_nodes.dat"
export NEKRS_CACHE_DIR="$PWD/.cache_gpu"
export OCCA_CACHE_DIR="$NEKRS_CACHE_DIR/occa"
export NEKRS_JITC_NTHREADS=1
mkdir -p "$OCCA_CACHE_DIR"
: "${HPC_PMIX:?HiPerGator MPI module must set HPC_PMIX}"
unset SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE SLURM_CPU_BIND_VERBOSE
echo "[RE_FULL_T100] Uout-S1__Pin-S1__Re0100_GL32x32 Re_D=100 P=6 U_inf=1.00 Job=$SLURM_JOB_ID GPUs=2 Nodes=1 t_start=$(date -Is)"
echo "[RE_FULL_T100] Uout-S1__Pin-S1__Re0100_GL32x32 mesh=$PWD/np3d.re2 NEKRS_HOME=$NEKRS_HOME"
srun --mpi="$HPC_PMIX" --ntasks=2 --ntasks-per-node=2 --cpus-per-task=8 \
   --gpus-per-task=1 --cpu-bind=none \
   "$NEKRS_HOME/bin/nekrs" --setup np3d.par --device-id 0
echo "[RE_FULL_T100] Done Uout-S1__Pin-S1__Re0100_GL32x32 $(date -Is)"
