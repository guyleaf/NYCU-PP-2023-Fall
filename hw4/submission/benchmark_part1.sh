#!/usr/bin/env bash
set -e

part1=("pi_block_linear" "pi_block_tree" "pi_gather" "pi_nonblock_linear" "pi_reduce" "pi_one_side")

./build.sh part1

if [[ -n ${1} ]]; then
	part1=("${1}")
fi

for i in "${part1[@]}"; do
	echo "${i}"
	echo "reference"
	mpirun --hostfile part1_hosts -np 4 "/home/HW4/ref/${i}" 1000000000 2>/dev/null | tail -n2
	echo "student"
	mpirun --hostfile part1_hosts -np 4 "./part1/${i}.out" 1000000000 2>/dev/null | tail -n2
	echo "=========================="
done
