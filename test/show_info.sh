#!/bin/bash

set -x

SCRIPT_DIR=$(cd $(dirname $0); pwd)
cd $SCRIPT_DIR

# show processor information
cat /proc/cpuinfo

# show memory information
free -m
