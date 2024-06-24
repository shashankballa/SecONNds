. scripts/common.sh

if [ $# -lt 1 ] || [ $# -gt 2 ]
then
    echo -e "${RED}Invalid number of arguments.${NC}"
    echo "Usage: run-all-fworks-client.sh [sqnet|resnet50|densenet121] [(optional) log-file tag]"
    exit 1
fi

if ! contains "sqnet resnet50 densenet121" $1; then
    echo -e "Usage: run-all-fworks-client.sh ${RED}[sqnet|resnet50|densenet121]${NC} [(optional) log-file tag]"
    exit 1
fi

bash scripts/run-client.sh SCI_HE $1 $2
bash scripts/run-client.sh seconnds_p $1 $2
bash scripts/run-client.sh cheetah $1 $2
bash scripts/run-client.sh seconnds_2 $1 $2
