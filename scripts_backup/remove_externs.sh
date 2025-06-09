echo "Starting removal of extern template declarations and macros from include/externs/"

# Define a function to process each file
process_file() {
    local file_path="$1"
    echo "Processing $file_path..."

    # Remove lines containing 'extern template'
    # Use a temporary file for sed in-place editing
    sed '/extern template/d' "$file_path" > "$file_path.tmp" && mv "$file_path.tmp" "$file_path"

    # Remove lines containing TFHEPP_EXPLICIT_INSTANTIATION macros
    # This regex matches lines that start with TFHEPP_EXPLICIT_INSTANTIATION, possibly with leading whitespace
    sed '/^[[:space:]]*TFHEPP_EXPLICIT_INSTANTIATION_/d' "$file_path" > "$file_path.tmp" && mv "$file_path.tmp" "$file_path"

    # Specific cleanup for include/externs/gate.hpp based on previous observations
    if [[ "$file_path" == "include/externs/gate.hpp" ]]; then
        sed '/INST(lvl1param);/d' "$file_path" > "$file_path.tmp" && mv "$file_path.tmp" "$file_path"
        sed '/INST(lvl0param);/d' "$file_path" > "$file_path.tmp" && mv "$file_path.tmp" "$file_path"
    fi

    # Specific cleanup for include/externs/trgsw.hpp based on previous observations
    # Remove lines like INST(lvl1param); that are not covered by the general TFHEPP_EXPLICIT_INSTANTIATION_ macro
    if [[ "$file_path" == "include/externs/trgsw.hpp" ]]; then
        # This is a bit more complex as the INST macro itself is defined and then undef'd multiple times.
        # The simple sed '/INST(lvl1param);/d' might be too greedy if INST is used differently elsewhere.
        # However, given the context of extern template removal, these specific instantiations are targeted.
        # A more robust way would be to parse, but for now, a targeted sed should work for this specific file structure.
        # Looking at the file, it seems these are direct calls after a #define INST(...) extern template ...
        # So, removing them when they appear alone on a line should be safe in this context.
        sed '/^INST(lvl1param);$/d' "$file_path" > "$file_path.tmp" && mv "$file_path.tmp" "$file_path"
        # also for lvl0param if it exists and follows the same pattern
        sed '/^INST(lvl0param);$/d' "$file_path" > "$file_path.tmp" && mv "$file_path.tmp" "$file_path"
    fi

    echo "Done processing $file_path."
}

# Export the function so it can be used by find -exec
export -f process_file

# Find all .hpp files in include/externs/ and process them
find include/externs/ -type f -name "*.hpp" -exec bash -c 'process_file "$0"' {} \;

echo "Removal of extern template declarations and macros from include/externs/ complete."

# Additionally, there might be files directly in include/ that use such macros,
# although 'externs' was the primary location.
# For example, if a file in 'include/' directly used 'extern template' or one of the TFHEPP_EXPLICIT_INSTANTIATION macros.
# A quick grep can check this.
echo "Searching for remaining extern template or TFHEPP_EXPLICIT_INSTANTIATION macros in include/ (excluding externs)..."
if grep -Er --include='*.hpp' --exclude-dir='externs' 'extern template|TFHEPP_EXPLICIT_INSTANTIATION_' include/; then
    echo "WARNING: Found remaining extern template or TFHEPP_EXPLICIT_INSTANTIATION macros in include/ (excluding externs). Manual review might be needed."
else
    echo "No remaining extern template or TFHEPP_EXPLICIT_INSTANTIATION macros found in include/ (excluding externs)."
fi
