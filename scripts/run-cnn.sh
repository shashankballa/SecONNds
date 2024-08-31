. scripts/common.sh

BASE_FWORK=cheetah # Base Framework
NTHREADS=$NUM_THREADS
SNN=0
NTT=0
NTRIPS=0
ROLE=client
ROLENUM=2
PRIVINP=0

echo -e "=============================================================================="
echo -e " "

if [ $# -lt 3 ]; then
  echo -e "${RED}Invalid number of arguments.${NC}"
  echo -e "Usage: run-cnn.sh [server|client] [seconnds_2|seconnds_p|cheetah|SCI_HE] [sqnet|resnet50|densenet121] [(optional) log-file tag]"
  exit 1
fi

if ! contains "server client" $1; then
  echo -e "Usage: run-cnn.sh ${RED}[server|client]${NC} [seconnds_2|seconnds_p|cheetah|SCI_HE] [sqnet|resnet50|densenet121] [(optional) log-file tag]"
  exit 1
fi

if ! contains "cheetah SCI_HE seconnds_2 seconnds_p" $2; then
  echo -e "Usage: run-cnn.sh [server|client] ${RED}[seconnds_2|seconnds_p|cheetah|SCI_HE]${NC} [sqnet|resnet50|densenet121] [(optional) log-file tag]"
  exit 1
fi

if ! contains "sqnet resnet50 densenet121" $3; then
  echo -e "Usage: run-cnn.sh [server|client] [seconnds_2|seconnds_p|cheetah|SCI_HE] ${RED}[sqnet|resnet50|densenet121]${NC} [(optional) log-file tag]"
  exit 1
fi

ROLE=$1

if [ "$ROLE" = "server" ]; then
  ROLENUM=1
  PRIVINP=pretrained/$3_model_scale12.inp
else
  ROLENUM=2
  PRIVINP=pretrained/$3_input_scale12_pred*.inp
fi

if [ "$2" = "SCI_HE" ]; then
  BASE_FWORK=SCI_HE
fi

if [ "$2" = "seconnds_2" ]; then
  BASE_FWORK=cheetah
  SNN=1
  NTT=1
  if [ "$3" = "sqnet" ]; then
    NTRIPS=$NTRIPS_2_SQNET
  fi
fi

if [ "$2" = "seconnds_p" ]; then
  BASE_FWORK=SCI_HE
  SNN=1
  NTT=1
  if [ "$3" = "sqnet" ]; then
    NTRIPS=$NTRIPS_P_SQNET
  fi
fi

# create a data/ to store the Ferret output
mkdir -p data

echo -e "$ROLE: Running ${GREEN}$3${NC} with ${GREEN}$2${NC}..."
echo -e " "

if [[ "$*" == *"--debug"* ]]; then
  gdb build/bin/$3-$BASE_FWORK -ex "run r=$ROLENUM k=$FXP_SCALE ell=$SS_BITLEN nt=$NTHREADS ip=$SERVER_IP p=$SERVER_PORT snn=$SNN ntrips=$NTRIPS csize=$CSIZE < $PRIVINP"
else
  if [ "$4" = "no_log" ]; then
    cat $PRIVINP | build/bin/$3-$BASE_FWORK r=$ROLENUM k=$FXP_SCALE ell=$SS_BITLEN nt=$NTHREADS ip=$SERVER_IP p=$SERVER_PORT snn=$SNN ntrips=$NTRIPS csize=$CSIZE
  else
    mkdir -p $LOGS_DIR
    if [ -z "$4" ]; then
      LOGFILE="$LOGS_DIR/$3-$2-$ROLE.log"
    else
      LOGFILE="$LOGS_DIR/$4-$3-$2-$ROLE.log"
    fi
    echo -e "Date: $(date)" > $LOGFILE
    echo -e "$ROLE: Running $3 with $2..." >> $LOGFILE
    echo -e " " >> $LOGFILE
    cat $PRIVINP | build/bin/$3-$BASE_FWORK r=$ROLENUM k=$FXP_SCALE ell=$SS_BITLEN nt=$NTHREADS ip=$SERVER_IP p=$SERVER_PORT snn=$SNN ntrips=$NTRIPS csize=$CSIZE >> $LOGFILE
    echo -e "Done! Log saved to:"
    echo -e "${GREEN}$LOGFILE${NC}"
  fi
fi

echo -e " "
echo -e "=============================================================================="
