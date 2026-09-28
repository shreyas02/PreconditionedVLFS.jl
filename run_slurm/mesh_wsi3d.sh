#!/bin/bash
# Run this script from the PreconditionedVLFS.jl directory

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

source ./run_slurm/env.sh

GEO_FILE="${GEO_FILE:-scripts/mesh_wsi_3d_2.geo}"
MESH_FILE="${MESH_FILE:-scripts/mesh_wsi_3d_2.msh}"
GMSH_THREADS="${GMSH_THREADS:-16}"
GMSH_PARTITIONS="${GMSH_PARTITIONS:-64}"

mkdir -p "${ROOT_DIR}/slurm_jobs"

echo "Meshing ${GEO_FILE} -> ${MESH_FILE}"
echo "Gmsh threads: ${GMSH_THREADS}; partitions: ${GMSH_PARTITIONS}"

srun --ntasks=1 --cpus-per-task="${GMSH_THREADS}" \
    gmsh "${GEO_FILE}" -3 -nt "${GMSH_THREADS}" -part "${GMSH_PARTITIONS}" -o "${MESH_FILE}" \
    > "${ROOT_DIR}/slurm_jobs/mesh_wsi_3d_2.log" 2>&1

echo "WSI 3D mesh completed."
