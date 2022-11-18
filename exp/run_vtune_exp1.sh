#!/bin/bash

# Run experiments that place all data structures on node 0 DRAM vs node 1 DRAM.

INPUT_DIR="/ssd1/songxin8/thesis/genomics/input-datasets/bwa-mem2-ert/"
NUM_THREADS=16
RESULT_DIR="/ssd1/songxin8/thesis/genomics/vtune/BWA-MEM2-ert/exp1/"


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

  VTUNE_MEMACC_COMMON="/opt/intel/oneapi/vtune/latest/bin64/vtune -collect memory-access \
      -knob sampling-interval=10 -knob analyze-mem-objects=true -knob analyze-openmp=true \
      -data-limit=10000 -result-dir ${RESULT_DIR}/${OUTFILE}_memacc"

  VTUNE_HOTSPOT_COMMON="/opt/intel/oneapi/vtune/latest/bin64/vtune -collect hotspots \
      -data-limit=10000 -result-dir ${RESULT_DIR}/${OUTFILE}_hotspot"

  VTUNE_UARCH_COMMON="/opt/intel/oneapi/vtune/latest/bin64/vtune -collect uarch-exploration \
      -knob sampling-interval=10 -knob collect-memory-bandwidth=true -data-limit=10000 -result-dir ${RESULT_DIR}/${OUTFILE}_uarch"

  # first 200M reads
  READS_FILE="${INPUT_DIR}/SRR7733443_200m_1.fastq"
  NUMACTL_COMMON="/usr/bin/numactl --membind=${MEMBIND} --cpunodebind=${CPUNODEBIND}"

  # Memory access analysis
  echo "Running memacc analysis. Log is at ${RESULT_DIR}/${OUTFILE}_memacc_log"
  ${VTUNE_MEMACC_COMMON} -- ${NUMACTL_COMMON} \
    ../bwa-mem2 mem -Y -K 100000000 -t ${NUM_THREADS} -Z ${INPUT_DIR}/human_g1k_v37_ert \
    $READS_FILE -o ${OUTFILE}.sam &> ${RESULT_DIR}/${OUTFILE}_memacc_log
  
  clean_cache

  # Hotspot analysis
  echo "Running hotspot analysis. Log is at ${RESULT_DIR}/${OUTFILE}_hotspot_log"
  ${VTUNE_HOTSPOT_COMMON} -- ${NUMACTL_COMMON} \
    ../bwa-mem2 mem -Y -K 100000000 -t ${NUM_THREADS} -Z ${INPUT_DIR}/human_g1k_v37_ert \
    $READS_FILE -o ${OUTFILE}.sam &> ${RESULT_DIR}/${OUTFILE}_hotspot_log

  clean_cache

  # uarch analysis
  echo "Running uarch analysis. Log is at ${RESULT_DIR}/${OUTFILE}_uarch_log"
  ${VTUNE_UARCH_COMMON} -- ${NUMACTL_COMMON} \
    ../bwa-mem2 mem -Y -K 100000000 -t ${NUM_THREADS} -Z ${INPUT_DIR}/human_g1k_v37_ert \
    $READS_FILE -o ${OUTFILE}.sam &> ${RESULT_DIR}/${OUTFILE}_uarch_log

  clean_cache
}


##############
# Script start
##############
trap clean_up SIGHUP SIGINT SIGTERM

[[ $EUID -ne 0 ]] && echo "This script must be run using sudo or as root." && exit 1

mkdir -p $RESULT_DIR

disable_autonuma 
clean_cache
# For all local, place all processes and data on node 1 memory
run_SRR7733443_1 "SRR7733443_1_${NUM_THREADS}threads_alllocal" "ALL_LOCAL"

disable_autonuma 
clean_cache                                                                
run_SRR7733443_1 "SRR7733443_1_${NUM_THREADS}threads_allremote" "ALL_REMOTE"

