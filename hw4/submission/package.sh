#!/usr/bin/env bash

id=$(id -u -n)
zip -x "*/.*" -r "HW4_${id}.zip" part1 part2
