#!/bin/bash

set -u

##############
# NUMA config
#############
cpunode="0"
memnode="0"

##############
# parameter for the bench
##############
logN=15
MM=(2 7 15)
#MM=(15)
hh=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15)
#hh=(1)
L=9
ll=(1 4 9)
kk=(1 5 10)
delta=40 # scaling factor
#d= # maximal number of digits for key-switching (an integer that divides L+1)
dd=(1 2 10)
T=100

GITREV=$(git rev-parse --short HEAD)


function usage() {
	echo "Usage: $0 [--thread] <normal|hexl>"
}

function eval_cmd() {
	eval numactl --cpunodebind 1 --membind 1 $@
}

function run_test() {
	mode=$1
	shift 1
	for M in "${MM[@]}"; do
		for h in `seq 1 $M`; do
			for l in "${ll[@]}"; do
				for d in "${dd[@]}"; do
					echo "[+] $logN $M $h $L $l $delta $d $T"
					if [ $M = 15 ] && [ $h = 1 ]; then
						echo "[+] Skip M=15, h=1 case due to OOM"
						continue
					fi
					eval_cmd "build-$mode/src/bench_unroll $logN $M $h $L $l $delta $d $T $*"
				done
			done
		done
	done
}

# main <mode> <time> <any flags that passed to run_test command>
function main() {
	mode=$1
	if [ $# -lt 2 ]; then
		time=$(date "+%Y%m%d_%H%M%S")
	else
		time=$2
	fi
	shift 2
	case $mode in
		normal)
			rm -rf /tmp/result/Trace{Runtime,Memory}
			mkdir -p /tmp/result/Trace{Runtime,Memory}
			run_test normal $*
			mv /tmp/result/TraceRuntime /tmp/result/TraceRuntime_${GITREV}_normal_$time
			mv /tmp/result/TraceMemory /tmp/result/TraceMemory_${GITREV}_normal_$time
			;;
		hexl)
			rm -rf /tmp/result/Trace{Runtime,Memory}
			mkdir -p /tmp/result/Trace{Runtime,Memory}
			run_test hexl $*
			mv /tmp/result/TraceRuntime /tmp/result/TraceRuntime_${GITREV}_hexl_$time
			mv /tmp/result/TraceMemory /tmp/result/TraceMemory_${GITREV}_hexl_$time
			;;
		*)
			;;
	esac
}

if [ $# -lt 1 ]; then
	usage
	exit -1
fi

time=$(date "+%Y%m%d_%H%M%S")

flag=""

if [ $1 = "--thread" ]; then
	flag="--thread"
fi

while [ $# -gt 0 ]
do
	echo "Run $1"
	main $1 $time $flag
	shift
done
