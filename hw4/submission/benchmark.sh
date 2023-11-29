#!/usr/bin/env bash
set -e

part1=("pi_block_linear" "pi_block_tree" "pi_gather" "pi_nonblock_linear" "pi_reduce" "pi_one_side")

for i in "${part1[@]}"; do
	echo "${i}"
	echo "reference"
	mpirun --hostfile hosts -np 4 "/home/HW4/ref/${i}" 1000000000
	echo "student"
	mpirun --hostfile hosts -np 4 "./part1/${i}.out" 1000000000
	echo "=========================="
done