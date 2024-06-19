. scripts/common.sh

SNN=0
# NTT=0

if [ $# -lt 2 ] || [ $# -gt 3 ]
then
  echo -e "${RED}Invalid number of arguments.${NC}"
  echo "Usage: run-client.sh [cheetah|seconnds] [sqnet|resnet50|densenet121] [(optional) log-file tag]"
exit 1
fi

if ! contains "cheetah seconnds" $1; then
  echo -e "Usage: run-client.sh ${RED}[cheetah|seconnds]${NC} [sqnet|resnet50|densenet121] [(optional) log-file tag]"
exit 1
fi

if ! contains "sqnet resnet50 densenet121" $2; then
  echo -e "Usage: run-client.sh [cheetah|seconnds] ${RED}[sqnet|resnet50|densenet121]${NC} [(optional) log-file tag]"
exit 1
fi

# Determine the log filename
if [ -z "$3" ]; then
  LOGFILE="$LOGS_DIR/$2-$1_client.log"
else
  LOGFILE="$LOGS_DIR/$3-$2-$1_client.log"
fi

if [ "$1" = "seconnds" ]; then
  SNN=1
  # NTT=0
fi

# create a data/ to store the Ferret output
mkdir -p data
mkdir -p logs

echo -e "Running ${GREEN}build/bin/$2-cheetah${NC}, which might take a while...."
cat pretrained/$2_input_scale12_pred*.inp | build/bin/$2-cheetah r=2 k=$FXP_SCALE ell=$SS_BITLEN nt=$NUM_THREADS ip=$SERVER_IP p=$SERVER_PORT snn=$SNN > $LOGFILE
echo -e "Computation done, check out the log file ${GREEN}$LOGFILE${NC}"
