echo "Starting replacement of #include with #import..."

find_files() {
    find include/ src/ -type f \( -name "*.hpp" -o -name "*.cpp" \)
}

# Process each file
while IFS= read -r file_path; do
    echo "Processing $file_path..."
    # Use a temporary file for sed in-place editing to handle errors robustly
    tmp_file="$file_path.tmp"

    # Replace #include <header> with #import <header>
    # The regex captures the header name including any paths within the angle brackets
    sed -E 's/^([[:space:]]*)#include[[:space:]]*<([^>]+)>/\1#import <\2>/' "$file_path" > "$tmp_file" && mv "$tmp_file" "$file_path"

    # Replace #include "header" with #import "header"
    # The regex captures the header name including any paths within the double quotes
    sed -E 's/^([[:space:]]*)#include[[:space:]]*"([^"]+)"/\1#import "\2"/' "$file_path" > "$tmp_file" && mv "$tmp_file" "$file_path"

    # Note: #pragma once is intentionally left untouched.
    # Semicolons are no longer added by these sed commands.

    echo "Done processing $file_path."
done < <(find_files)

echo "Replacement of #include with #import complete."

# Verification (optional, but good for a spot check)
echo "Spot-checking a few files for #import statements..."
find_files | head -n 5 | while IFS= read -r file_path; do
    echo "--- Content of $file_path (first 10 lines) ---"
    head -n 10 "$file_path"
    echo "-------------------------------------------"
done
