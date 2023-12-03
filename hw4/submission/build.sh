#!/usr/bin/env bash
set -e

make clean
make
parallel-scp -h part1_hosts -r ~/HW4 ~
