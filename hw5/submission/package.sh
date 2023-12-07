#!/usr/bin/env bash

id=$(id -u -n)
zip "HW5_${id}.zip" kernel*.cu
