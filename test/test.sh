#!/bin/bash

set -eu

# thanks to https://qiita.com/Hayao0819/items/0e04b39b0804a0d16020
contains() {
    local i
    for i in $2{[@]}; do
        if [[ ${i} = ${1} ]]; then
            echo "TRUE"
            return
        fi
    done
    echo "FALSE"
    return
}

SCRIPT_DIR=$(cd $(dirname $0); pwd)
cd $SCRIPT_DIR

./show_info.sh

TEST_BINARIES=`find . -maxdepth 1 -perm -111 -type f`
IGNORED_BINARIES=('')
for TEST_BINARY in $TEST_BINARIES
do
    RES=`contains $TEST_BINARY "${IGNORED_BINARIES[*]}"`
    if [[ $RES == "TRUE" ]]; then
        continue
    fi
    echo "==== Testing $TEST_BINARY ===="
    ./$TEST_BINARY
done
