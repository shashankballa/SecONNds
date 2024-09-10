. scripts/common.sh

BUILD_MODE=Release
TRIP_TRIALS=0

for deps in eigen3 emp-ot emp-tool hexl SEAL-4.0
do
  if [ ! -d $BUILD_DIR/include/$deps ] 
  then
	echo -e "${RED}$deps${NC} seems absent in ${BUILD_DIR}/include/, please re-run scripts/build-deps.sh"
	exit 1
  fi
done

for deps in zstd.h 
do
  if [ ! -f $BUILD_DIR/include/$deps ] 
  then
	echo -e "${RED}$deps${NC} seems absent in ${BUILD_DIR}/include/, please re-run scripts/build-deps.sh"
	exit 1
  fi
done

cd $BUILD_DIR/

if [[ $* == *"--debug"* ]]; then
  BUILD_MODE=Debug
fi

if [[ $* == *"--trip_trials"* ]]; then
  TRIP_TRIALS=1
fi

cmake .. -DCMAKE_BUILD_TYPE=$BUILD_MODE -DSCI_BUILD_NETWORKS=ON -DSCI_BUILD_TESTS=ON \
          -DOPENSSL_ROOT_DIR=/usr/local/opt/openssl -DCMAKE_PREFIX_PATH=$BUILD_DIR -DUSE_APPROX_RESHARE=ON \
          -DRUN_TRIP_TRIALS=$TRIP_TRIALS

for net in resnet50 sqnet densenet121
do
    make ${net}-cheetah -j16
    make ${net}-SCI_HE -j16
done
