#!/usr/bin/env bash
set -e

processes=(2 4 6 8)
data=($(ls /home/.grade/HW4/data-set/data*))

./build.sh part2

for i in "${data[@]}"; do
	echo "${i}"
	for j in "${processes[@]}"; do
		echo "Run in ${j} processes"
		mpirun --hostfile part2_hosts -np "${j}" "./part2/matmul" <"${i}" 2>/dev/null | tail -n"${j}"
		echo "=========================="
		sleep 1
	done
done
