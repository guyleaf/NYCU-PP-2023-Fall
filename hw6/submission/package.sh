#!/usr/bin/env bash

id=$(id -u -n)
zip "HW6_${id}.zip" hostFE.c kernel.cl url.txt
