. scripts/common.sh

cd $DEPS_DIR/troy-nova
bash scripts/test_conv2d.sh -ntt -logdir="$LOGS_DIR" $*
bash scripts/test_conv2d.sh -logdir="$LOGS_DIR" $*