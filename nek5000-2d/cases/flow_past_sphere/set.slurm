#!/bin/bash
#SBATCH --job-name=U14G048
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%j.log
#SBATCH --error=slurm_%j.err

set -euo pipefail
fatal(){ echo "FATAL: $*" >&2; exit 2; }

if command -v module >/dev/null 2>&1; then
  module purge || true
  module load ufrc || true
  module load gcc/14.2 2>/dev/null || module load gcc/14.2.0
  module load openmpi/5.0.10
fi

export OMPI_MCA_coll_hcoll_enable=0
export FC="mpifort -mcmodel=large"
export LD="mpifort -mcmodel=large"
export CC="mpicc"
export FFLAGS="-fallow-argument-mismatch -std=legacy -w"

case_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"
case_dir="$(cd "$case_dir" && pwd)"
pkg_dir="$(cd "$case_dir/.." && pwd)"
cd "$case_dir"

bash ./prepare_mesh.sh

NEK_BASE="${NEK_BASE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
PRIVATE_NEK="$case_dir/.nek_uniformmesh"
rm -rf "$PRIVATE_NEK"
mkdir -p "$PRIVATE_NEK"
for item in "$NEK_BASE"/*; do
  b="$(basename "$item")"
  [[ "$b" == "run" || "$b" == "core" ]] && continue
  ln -s "$item" "$PRIVATE_NEK/$b"
done
cp -a "$pkg_dir/ACTIVE_CORE_BK/core" "$PRIVATE_NEK/core"
export NEK_SOURCE_ROOT="$PRIVATE_NEK"

rm -rf obj nek5000 set.f makefile .makefile .state compiler.out build.log
rm -f SESSION.NAME set0.f????? set.nek5000
rm -f drag.dat stage35AD_cv_drag.dat stage35AR_forces.dat
rm -f stage35AD_signature.dat stage35AD_coeff8.dat
rm -f stage23_mapping_audit.dat stage23_operator_audit.dat
rm -f stage24_dynamic_diag.dat stage27_meanfree_diag.dat
rm -f stage27_postsolve_diag.dat stage35BQ_lift_diag.dat
mkdir -p obj

printf 'set\n%s/\n' "$case_dir" > SESSION.NAME
echo N | bash ./makenek set
[[ -x ./nek5000 ]] || fatal "build failed"

if srun --mpi=list 2>&1 | grep -q pmix_v4; then
  srun --mpi=pmix_v4 -n 8 ./nek5000
elif srun --mpi=list 2>&1 | grep -q pmix; then
  srun --mpi=pmix -n 8 ./nek5000
else
  mpiexec -np 8 ./nek5000
fi
