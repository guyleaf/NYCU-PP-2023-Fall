#!/usr/bin/env bash
set -e

make clean
make "${@:1}"
parallel-scp -h part1_hosts -r ~/HW4 ~
