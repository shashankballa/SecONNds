. scripts/common.sh

BASE_FWORK=cheetah # Base Framework
NTHREADS=$NUM_THREADS
SNN=0
# NTT=0
NTRIPS=0

if [ $# -lt 2 ] || [ $# -gt 3 ]
then
  echo -e "${RED}Invalid number of arguments.${NC}"
  echo "Usage: run-client.sh [cheetah|SCI_HE|seconnds] [sqnet|resnet50|densenet121] [(optional) log-file tag]"
exit 1
fi

if [ "$1" = "SCI_HE" ]; then
  BASE_FWORK=SCI_HE
fi

if [ "$1" = "seconnds_2" ]; then
  BASE_FWORK=cheetah
  SNN=1
  # NTT=0
  if [ "$2" = "sqnet" ]; then
    NTRIPS=$NTRIPS_2_SQNET
  fi
fi

if [ "$1" = "seconnds_p" ]; then
  BASE_FWORK=SCI_HE
  SNN=1
  # NTT=0
  if [ "$2" = "sqnet" ]; then
    NTRIPS=$NTRIPS_P_SQNET
  fi
fi

if ! contains "cheetah SCI_HE seconnds_2 seconnds_p" $1; then
  echo -e "Usage: run-client.sh ${RED}[cheetah|seconnds]${NC} [sqnet|resnet50|densenet121] [(optional) log-file tag]"
exit 1
fi

if ! contains "sqnet resnet50 densenet121" $2; then
  echo -e "Usage: run-client.sh [cheetah|seconnds] ${RED}[sqnet|resnet50|densenet121]${NC} [(optional) log-file tag]"
exit 1
fi

# create a data/ to store the Ferret output
mkdir -p data


echo -e "CLIENT: Running ${GREEN}$2-$BASE_FWORK${NC} with ${GREEN}$1${NC}..."
echo -e " "

if [ "$3" = "no_log" ]; then
  cat pretrained/$2_input_scale12_pred*.inp | build/bin/$2-$BASE_FWORK r=2 k=$FXP_SCALE ell=$SS_BITLEN nt=$NTHREADS ip=$SERVER_IP p=$SERVER_PORT snn=$SNN ntrips=$NTRIPS csize=$CSIZE
else
  mkdir -p $LOGS_DIR
  if [ -z "$3" ]; then
    LOGFILE="$LOGS_DIR/$2-$1-client.log"
  else
    LOGFILE="$LOGS_DIR/$3-$2-$1-client.log"
  fi
  echo -e "Date: $(date)" > $LOGFILE
  echo -e "CLIENT: Running $2 with $1..." >> $LOGFILE
  echo -e " " >> $LOGFILE
  cat pretrained/$2_input_scale12_pred*.inp | build/bin/$2-$BASE_FWORK r=2 k=$FXP_SCALE ell=$SS_BITLEN nt=$NTHREADS ip=$SERVER_IP p=$SERVER_PORT snn=$SNN ntrips=$NTRIPS csize=$CSIZE >> $LOGFILE
  echo -e "Done! Log saved to ${GREEN}$LOGFILE${NC}"
fi