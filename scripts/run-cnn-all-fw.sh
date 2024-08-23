. scripts/common.sh

if [ $# -lt 2 ] || [ $# -gt 3 ]
then
    echo -e "${RED}Invalid number of arguments.${NC}"
    echo "Usage: run-cnn-all-fw [server|client] [sqnet|resnet50|densenet121] [(optional) log-file tag]"
    exit 1
fi

if ! contains "server client" $1; then
    echo -e "Usage: run-cnn-all-fw ${RED}[server|client]${NC} [sqnet|resnet50|densenet121] [(optional) log-file tag]"
    exit 1
fi

if ! contains "sqnet resnet50 densenet121" $2; then
    echo -e "Usage: run-cnn-all-fw [server|client] ${RED}[sqnet|resnet50|densenet121]${NC} [(optional) log-file tag]"
    exit 1
fi

for FW in SCI_HE cheetah seconnds_p seconnds_2
do
    bash scripts/run-cnn.sh $1 $FW $2 $3
done