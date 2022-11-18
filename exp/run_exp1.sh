#!/bin/bash

# Run experiments that place all data structures on node 0 DRAM vs node 1 DRAM.

INPUT_DIR="/ssd1/songxin8/thesis/genomics/input-datasets/bwa-mem2-ert/"
NUM_THREADS=16
RESULT_DIR="exp1"

clean_up () {
    echo "Cleaning up. Kernel PID is $EXE_PID, numastat PID is $NUMASTAT_PID."
    # Perform program exit housekeeping
    kill $NUMASTAT_PID
    kill $EXE_PID
    kill $TOP_PID
    exit
}

clean_cache () { 
  echo "Clearing caches..."
  # clean CPU caches
  ./tools/clear_cpu_cache
  # clean page cache
  echo 3 > /proc/sys/vm/drop_caches
}

disable_autonuma () {
  # turn off both numa
  sudo service numad stop
  NUMAD_OUT=$(systemctl is-active numad)
  echo "numad service is now $NUMAD_OUT (should be not active)"

  echo 0 > /proc/sys/vm/zone_reclaim_mode
  echo 0 > /proc/sys/kernel/numa_balancing
  NUMA_BALANCING=$(cat /proc/sys/kernel/numa_balancing)
  echo "numa_balancing is now $NUMA_BALANCING (should be 0)"
}

run_SRR7733443_1() { 
  OUTFILE=$1 #first argument
  CONFIG=$2

  if [[ "$CONFIG" == "ALL_LOCAL" ]]; then
    # All local config: place both data and compute on node 1
    CPUNODEBIND=1
    MEMBIND=1
  elif [[ "$CONFIG" == "ALL_REMOTE" ]]; then
    # Emulating CXL
    CPUNODEBIND=0
    MEMBIND=1
  else
    echo "Error! Undefined configuration $CONFIG"    
    exit 1
  fi

  echo "Start" > $OUTFILE
  
  # first 200M reads
  READS_FILE="${INPUT_DIR}/SRR7733443_200m_1.fastq"

  echo "/usr/bin/time -v /usr/bin/numactl --membind=${MEMBIND} --cpunodebind=${CPUNODEBIND} \
    ./bwa-mem2 mem -Y -K 100000000 -t ${NUM_THREADS} -Z ${INPUT_DIR}/human_g1k_v37_ert \
    $READS_FILE -o ${OUTFILE}.sam &>> $OUTFILE &"

  /usr/bin/time -v /usr/bin/numactl --membind=${MEMBIND} --cpunodebind=${CPUNODEBIND} \
    ../bwa-mem2 mem -Y -K 100000000 -t ${NUM_THREADS} -Z ${INPUT_DIR}/human_g1k_v37_ert \
    $READS_FILE -o ${OUTFILE}.sam &>> $OUTFILE &

  TIME_PID=$! 
  EXE_PID=$(pgrep -P $TIME_PID)

  echo "EXE PID is ${EXE_PID}"
  echo "start" > ${OUTFILE}_numastat
  while true; do numastat -p $EXE_PID >> ${OUTFILE}_numastat; sleep 5; done &
  NUMASTAT_PID=$!
  top -b -d 10 -1 -p $EXE_PID > ${OUTFILE}_top_log &
  TOP_PID=$!

  echo "Waiting for kernel to complete (PID is ${EXE_PID}). numastat is logged into ${OUTFILE}_numastat, PID is ${NUMASTAT_PID}. Top log PID is ${TOP_PID}" 
  wait $TIME_PID
  echo "GAP kernel complete."
  kill $NUMASTAT_PID
  kill $TOP_PID
}


##############
# Script start
##############
trap clean_up SIGHUP SIGINT SIGTERM

[[ $EUID -ne 0 ]] && echo "This script must be run using sudo or as root." && exit 1

mkdir -p $RESULT_DIR

#disable_autonuma 
#clean_cache
## For all local, place all processes and data on node 1 memory
#run_SRR7733443_1 "${RESULT_DIR}/SRR7733443_1_${NUM_THREADS}threads_alllocal" "ALL_LOCAL"

disable_autonuma 
clean_cache                                                                
# To emulate CXL, place processes on node 0
run_SRR7733443_1 "${RESULT_DIR}/SRR7733443_1_${NUM_THREADS}threads_allremote" "ALL_REMOTE"

