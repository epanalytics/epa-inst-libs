#!/usr/bin/env bash

pebildir="/home/ewoodruff/bitbucket/pebil"
export PEBIL_ROOT=$pebildir
echo "***** initializing environment using PEBIL_ROOT=$PEBIL_ROOT"
export LD_LIBRARY_PATH=${PEBIL_ROOT}/lib:${PEBIL_ROOT}/external/ReuseDistance:${LD_LIBRARY_PATH}
export PATH=${PEBIL_ROOT}/bin:${PEBIL_ROOT}/scripts:${PEBIL_ROOT}/external/udis86-1.7/udcli:${PATH}

