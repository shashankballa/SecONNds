. scripts/common.sh

SNN=0
NTT=0
BASE_FWORK=cheetah # Base Framework
NTHREADS=$NUM_THREADS

if [ $# -lt 2 ] || [ $# -gt 3 ]
then
  echo -e "${RED}Invalid number of arguments.${NC}"
  echo "Usage: run-client.sh [cheetah|SCI_HE|seconnds] [sqnet|resnet50|densenet121] [(optional) log-file tag]"
exit 1
fi

# if the first argument is "SCI_HE" then set base Framework to "SCI_HE"
if [ "$1" = "SCI_HE" ]; then
  BASE_FWORK=SCI_HE
  NTHREADS=4
fi

if [ "$1" = "seconnds" ]; then
  SNN=1
  NTT=1
fi

if ! contains "cheetah SCI_HE seconnds" $1; then
  echo -e "Usage: run-client.sh ${RED}[cheetah|seconnds]${NC} [sqnet|resnet50|densenet121] [(optional) log-file tag]"
exit 1
fi

if ! contains "sqnet resnet50 densenet121" $2; then
  echo -e "Usage: run-client.sh [cheetah|seconnds] ${RED}[sqnet|resnet50|densenet121]${NC} [(optional) log-file tag]"
exit 1
fi

# Determine the log filename
if [ -z "$3" ]; then
  LOGFILE="$LOGS_DIR/$2-$1_server.log"
else
  LOGFILE="$LOGS_DIR/$3-$2-$1_server.log"
fi

# create a data/ to store the Ferret output
mkdir -p data
mkdir -p logs
# ls -lh pretrained/$2_model_scale12.inp
echo -e "SERVER: Running ${GREEN}$2${NC} with ${GREEN}$1${NC}..."
cat pretrained/$2_model_scale12.inp | build/bin/$2-$BASE_FWORK r=1 k=$FXP_SCALE ell=$SS_BITLEN nt=$NTHREADS p=$SERVER_PORT snn=$SNN ntt=$NTT > $LOGFILE
echo -e "Done! Log saved to ${GREEN}$LOGFILE${NC}"
