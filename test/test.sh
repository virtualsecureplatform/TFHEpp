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

if [[ $# -gt 0 ]]; then
    TEST_DIRS=()
    for TEST_DIR in "$@"
    do
        case "$TEST_DIR" in
            ./*|/*) TEST_DIRS+=("$TEST_DIR") ;;
            *) TEST_DIRS+=("./$TEST_DIR") ;;
        esac
    done
else
    TEST_DIRS=(
        './common'
        './tfhe'
        './bfv'
        './clpx'
        './ckks'
    )
fi
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
if [[ ${#TEST_BINARIES[@]} -eq 0 ]]; then
    echo "No test binaries found for: ${TEST_DIRS[*]}" >&2
    exit 1
fi
IGNORED_BINARIES=(
    # Argument-driven, long-running benchmark; run through its dedicated runner.
    './clpx/scheme_switch_comparison'
    './tfhe/largelut_cmux_benchmark'
    './tfhe/largelut_random_key_diagnostic'
    # Retained-checkpoint tools require external fixtures and explicit paths.
    './ckks/ckks_auto_key_reduce_exact'
    './ckks/ckks_narrow_auto_primary_projection'
    './ckks/ckks_narrow_auxiliary_functions'
    './ckks/ckks_narrow_auxiliary_probe'
    './ckks/ckks_narrow_auxiliary_replay'
    './ckks/ckks_narrow_auxiliary_square'
    './ckks/ckks_narrow_basis_batch_checkpoint'
    './ckks/ckks_narrow_basis_factor_cache_checkpoint'
    './ckks/ckks_narrow_basis_replay'
    './ckks/ckks_narrow_dd_ntt'
    './ckks/ckks_narrow_dd_ntt_tree'
    './ckks/ckks_narrow_ll14_basis_ntt_fbs'
    './ckks/ckks_narrow_ll14_basis_setup'
    './ckks/ckks_narrow_ll14_lean37_fbs'
    './ckks/ckks_narrow_ll14_mixed_fbs'
    './ckks/ckks_narrow_ll14_native_auto_setup'
    './ckks/ckks_narrow_ll14_ntt_complex_fbs'
    './ckks/ckks_narrow_ll14_ntt_fbs'
    './ckks/ckks_narrow_ll14_packed_complex_fbs'
    './ckks/ckks_narrow_ll14_packed_fbs'
    './ckks/ckks_narrow_ll14_retained_fbs'
    './ckks/ckks_narrow_ll14_retuned_fbs'
    './ckks/ckks_narrow_ll14_retuned_native_setup'
    './ckks/ckks_narrow_ll14_retuned_setup'
    './ckks/ckks_narrow_ll14_two_layer_fbs'
    './ckks/ckks_narrow_packed_archive'
    './ckks/ckks_narrow_packed_primary_projection'
    './ckks/ckks_narrow_packed_range'
    './ckks/ckks_narrow_packed_reconstruction_exact'
    './ckks/ckks_narrow_selector_ntt_file'
    './ckks/ckks_ram_setup'
    './ckks/ckks_ram_fixture_check'
    './ckks/ckks_narrow_ll14_retuned_ram'
)
for TEST_BINARY in "${TEST_BINARIES[@]}"
do
    if contains "$TEST_BINARY" "${IGNORED_BINARIES[@]}"; then
        continue
    fi
    echo "==== Testing $TEST_BINARY ===="
    "$TEST_BINARY"
done
