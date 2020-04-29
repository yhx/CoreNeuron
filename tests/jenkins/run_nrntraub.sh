#!/usr/bin/bash

set -e
source ${JENKINS_DIR:-.}/_env_setup.sh
source ${JENKINS_DIR:-.}/_test_utils.sh

set -x

CORENRN_TYPE="$1"
MPI_RANKS="$2"
THREADS="$3"

cd $WORKSPACE/nrntraub

# Allocate 1 node
N=1

if [ "${CORENRN_TYPE}" = "SoA" ]; then
  export OMP_NUM_THREADS=${THREADS}
  n=${MPI_RANKS} bb5_run ./${CORENRN_TYPE}/special-core --mpi -d coredat --voltage=1000
elif [ "${CORENRN_TYPE}" = "AoS" ]; then
  bb5_run ./${CORENRN_TYPE}/special-core --mpi -d coredat --voltage=1000
else
    echo "Not a valid CORENRN_TYPE"
    exit 1
fi

sort -n -k'1,1' -k2 < out.dat | awk 'NR==1 { print; next } { printf "%.3f\t%d\n", $1, $2 }' > out_cn_${CORENRN_TYPE}.sorted

sdiff -s out_cn_${CORENRN_TYPE}.sorted out_nrn.sorted
