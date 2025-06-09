echo "Starting removal of semicolons from #import lines..."

find_files() {
    find include/ src/ -type f \( -name "*.hpp" -o -name "*.cpp" \)
}

# Process each file
while IFS= read -r file_path; do
    echo "Processing $file_path for semicolon removal..."
    # Use a temporary file for sed in-place editing
    tmp_file="$file_path.tmp"

    # Remove trailing semicolons from #import lines
    sed -E 's/^(#import[[:space:]]*[<"][^>"]+[>"])[[:space:]]*;/\1/' "$file_path" > "$tmp_file" && mv "$tmp_file" "$file_path"

    echo "Done processing $file_path."
done < <(find_files)

echo "Semicolon removal from #import lines complete."

# Verification (optional, but good for a spot check)
echo "Spot-checking a few files for #import statements (should have no semicolons)..."
find_files | head -n 5 | while IFS= read -r file_path; do
    echo "--- Content of $file_path (first 10 lines) ---"
    head -n 10 "$file_path"
    echo "-------------------------------------------"
done
