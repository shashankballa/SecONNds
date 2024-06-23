# Rebuilds the project
# Author: Shashank Balla

# Usage: bash rebuild_full.sh #1 (optional: remove the previous build if 1)

PREFIX_DIR=/home/shashank/codes/seconnds_cnn/build

# remove the previous build if the first argument is 1
if [ "$1" == "1" ]; then
    rm -rf build
fi

# create a new build directory
mkdir build

# enter the build directory
cd build

# run cmake
cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX_DIR -DCMAKE_PREFIX_PATH=$PREFIX_DIR

# make the project
make -j8

# install the project
make install

# exit the build directory
cd ..