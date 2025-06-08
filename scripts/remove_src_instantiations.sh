echo "Starting removal of explicit template instantiations and macros from src/"

# Define a function to process each file
process_src_file() {
    local file_path="$1"
    echo "Processing $file_path..."

    # Remove lines containing 'template void' or 'template class' etc. followed by TFHEPP_EXPLICIT_INSTANTIATION macro
    # This is tricky because 'template <...>' is valid C++. We are looking for explicit instantiations.
    # Explicit instantiations look like 'template void func<...>(...);' or 'template class Class<...>;'
    # The relevant files use macros like TFHEPP_EXPLICIT_INSTANTIATION_... to define these.
    # So, the primary target is removing those macro calls.

    # Remove lines containing TFHEPP_EXPLICIT_INSTANTIATION macros
    # This regex matches lines that start with TFHEPP_EXPLICIT_INSTANTIATION, possibly with leading whitespace
    sed '/^[[:space:]]*TFHEPP_EXPLICIT_INSTANTIATION_/d' "$file_path" > "$file_path.tmp" && mv "$file_path.tmp" "$file_path"

    # Also, remove the macro definitions themselves if they are within the .cpp files (unlikely but possible)
    # For instance, lines like "#define INST(P) template void funcname<P>(...)"
    # However, the main goal is to remove the *invocations* of these macros.
    # The previous sed command should handle lines like "TFHEPP_EXPLICIT_INSTANTIATION_...(INST)"

    # A more direct approach for the template instantiations themselves:
    # Remove lines that look like "template void functionName<...>()" or "template returnType functionName<...>()"
    # This needs to be done carefully to avoid removing template definitions.
    # Given the structure observed (macros expanding to 'template void/class ...'), removing the macro calls is safer.

    # Let's ensure that any direct instantiations like "template void MyFunction<...>(...);" are also removed.
    # These are often wrapped by the INST macros that were just removed by the sed command above.
    # For example, if a file had:
    # #define INST(P) template void my_function<P>();
    # INST(type1);
    # The `sed` for TFHEPP_EXPLICIT_INSTANTIATION might not get it if `INST` is called directly.
    # The files show patterns like:
    # #define INST(...) template ...
    # TFHEPP_MACRO(INST)
    # #undef INST
    # So removing the TFHEPP_MACRO calls is key.

    # What if there are files that do not use the TFHEPP_EXPLICIT_INSTANTIATION macros but still have
    # "template class" or "template struct" or "template void function" instantiations?
    # The current approach of removing TFHEPP_EXPLICIT_INSTANTIATION macros is based on the observed pattern.
    # If direct explicit instantiations exist outside this pattern, they would need a different rule.
    # For now, assuming the macros are the primary mechanism.

    echo "Done processing $file_path."
}

# Export the function so it can be used by find -exec
export -f process_src_file

# Find all .cpp files in src/ and process them
find src/ -type f -name "*.cpp" -exec bash -c 'process_src_file "$0"' {} \;

echo "Removal of explicit template instantiations and macros from src/ complete."

# Verify no TFHEPP_EXPLICIT_INSTANTIATION macros remain in src files
echo "Verifying removal of TFHEPP_EXPLICIT_INSTANTIATION macros in src/..."
if grep -Er --include='*.cpp' 'TFHEPP_EXPLICIT_INSTANTIATION_' src/; then
    echo "WARNING: Found remaining TFHEPP_EXPLICIT_INSTANTIATION macros in src/. Manual review might be needed."
else
    echo "Successfully removed all TFHEPP_EXPLICIT_INSTANTIATION macro invocations from src/ files."
fi
