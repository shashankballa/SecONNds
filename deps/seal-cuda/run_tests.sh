# Change to build directory
cd build

# Run ctest
ctest

# Run timing tests
test/timetest_seal
test/timetest

# Change back to original directory
cd ..