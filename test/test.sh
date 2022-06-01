#!/bin/bash

set -eux

SCRIPT_DIR=$(cd $(dirname $0); pwd)
cd $SCRIPT_DIR

TEST_BINARIES=`find . -maxdepth 1 -perm -111 -type f`
for TEST_BINARY in $TEST_BINARIES
do
    echo "Testing $TEST_BINARY"
    ./$TEST_BINARY
done