#!/bin/bash

set -eu

# thanks to https://qiita.com/Hayao0819/items/0e04b39b0804a0d16020
contains() {
    local needle=$1
    shift
    local i
    for i in "$@"; do
        if [[ ${i} = ${needle} ]]; then
            return 0
        fi
    done
    return 1
}

SCRIPT_DIR=$(cd $(dirname $0); pwd)
cd $SCRIPT_DIR

./show_info.sh

TEST_DIRS=(
    './common'
    './tfhe'
    './bfv'
    './clpx'
    './ckks'
)
TEST_BINARIES=()
for TEST_DIR in "${TEST_DIRS[@]}"
do
    if [[ ! -d ${TEST_DIR} ]]; then
        continue
    fi
    while IFS= read -r TEST_BINARY
    do
        TEST_BINARIES+=("${TEST_BINARY}")
    done < <(find "${TEST_DIR}" -maxdepth 1 -perm -111 -type f | sort)
done
IGNORED_BINARIES=(
    './tfhe/largelut_cmux_benchmark'
    './tfhe/largelut_random_key_diagnostic'
)
for TEST_BINARY in "${TEST_BINARIES[@]}"
do
    if contains "$TEST_BINARY" "${IGNORED_BINARIES[@]}"; then
        continue
    fi
    echo "==== Testing $TEST_BINARY ===="
    "$TEST_BINARY"
done
