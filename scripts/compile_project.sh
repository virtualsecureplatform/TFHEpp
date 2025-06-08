echo "Attempting to compile the project..."

# Clean build: remove existing build directory
if [ -d "build" ]; then
    echo "Removing existing build directory..."
    rm -rf build
fi

# Configure with CMake
echo "Configuring with CMake..."
cmake -B build . -G Ninja -DENABLE_TEST=ON
# Check if CMake configuration was successful
if [ $? -ne 0 ]; then
    echo "CMake configuration failed."
    exit 1
fi
echo "CMake configuration successful."

# Build with Ninja
echo "Building with Ninja..."
cd build
ninja
# Check if Ninja build was successful
if [ $? -ne 0 ]; then
    echo "Ninja build failed."
    exit 1
fi
echo "Ninja build successful."

echo "Compilation attempt finished."
