. scripts/common.sh

SNN=0

if [ ! $# -eq 2 ]
then
  echo -e "${RED}Please specify the network to run.${NC}"
  echo "Usage: run-server.sh [cheetah|seconnds] [sqnet/resnet50]"
else
  if ! contains "cheetah seconnds" $1; then
    echo -e "Usage: run-server.sh ${RED}[cheetah|seconnds]${NC} [sqnet|resnet50|densenet121]"
	exit 1
  fi

  if ! contains "sqnet resnet50 densenet121" $2; then
    echo -e "Usage: run-server.sh [cheetah|seconnds] ${RED}[sqnet|resnet50|densenet121]${NC}"
	exit 1
  fi

  if [ "$1" = "seconnds" ]; then
    SNN=1
  fi

  # create a data/ to store the Ferret output
  mkdir -p data
  mkdir -p logs
  ls -lh pretrained/$2_model_scale12.inp
  echo -e "Running ${GREEN}build/bin/$2-cheetah${NC}, which might take a while...."
  cat pretrained/$2_model_scale12.inp | build/bin/$2-cheetah r=1 k=$FXP_SCALE ell=$SS_BITLEN nt=$NUM_THREADS p=$SERVER_PORT snn=$SNN > logs/$1-$2_server.log
  echo -e "Computation done, check out the log file ${GREEN}$1-$2_server.log${NC}"
fi
