#!/usr/bin/env bash
set -e

make clean
make
parallel-scp -h hosts -r ~/HW4 ~
