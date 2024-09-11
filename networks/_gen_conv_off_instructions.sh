
KEEP_INTERMEDIATE_FILES=0

if [[ $* == *"--keep-files"* ]]; then
    KEEP_INTERMEDIATE_FILES=1
    # remove --keep-files from the arguments
    set -- "${@/--keep-files/}"
fi

# if less than 1 arguments are provided
if [ $# -lt 1 ]; then
    echo -e "${RED}Invalid number of arguments.${NC}"
    echo "Usage: extract_tmp_pts.sh <file_name>"
    exit 1
fi

INP_FNAME=$1

# function to generate output filename
generate_out_fname() {
    local INP_FNAME=$1
    local OUT_TAG=$2
    local OUT_FNAME="${INP_FNAME%.*}_$OUT_TAG.${INP_FNAME##*.}"
    echo "$OUT_FNAME"
}

OUT_TAG=0
python3 _extract_pattern.py $INP_FNAME "Conv2DWrapper.*?;" -outtag="$OUT_TAG" --flatten
OUT0_FNAME=$(generate_out_fname "$INP_FNAME" "$OUT_TAG")

OUT_TAG=1
python3 _extract_pattern.py $OUT0_FNAME "Conv2DWrapper\([^)]*, ([^,]+), [^,]+\);" -outtag="$OUT_TAG"\
    -prefix="std::vector<std::vector<seal::Plaintext>> " -suffix="_pts;" -offset=6
OUT1_FNAME=$(generate_out_fname "$OUT0_FNAME" "$OUT_TAG")

mv $OUT1_FNAME $(generate_out_fname "$INP_FNAME" "tmp_pts")

OUT_TAG=2
python3 _extract_pattern.py $OUT0_FNAME "Conv2DWrapper\(([^,]+(?:, [^,]+)*?), [^,]+\);" -outtag="$OUT_TAG"\
    -offset=6 #-prefix="ConvOfflineCheetah(conv_ntt, " -suffix=");" 

# output filename is the input filename with _res_2 appended before the extension
OUT2_FNAME=$(generate_out_fname "$OUT0_FNAME" "$OUT_TAG")

if [ $KEEP_INTERMEDIATE_FILES -eq 0 ]; then
    rm $OUT0_FNAME
fi

OUT_TAG=3
python3 _drop_index.py $OUT2_FNAME 13 -outtag="$OUT_TAG" -offset=5
OUT3_FNAME=$(generate_out_fname "$OUT2_FNAME" "$OUT_TAG")

if [ $KEEP_INTERMEDIATE_FILES -eq 0 ]; then
    rm $OUT2_FNAME
fi

OUT_TAG=4
python3 _modify_and_insert.py $OUT3_FNAME 13 14 -outtag="$OUT_TAG" \
    -prefix="" -suffix="_pts" -offset=3
OUT4_FNAME=$(generate_out_fname "$OUT3_FNAME" "$OUT_TAG")

if [ $KEEP_INTERMEDIATE_FILES -eq 0 ]; then
    rm $OUT3_FNAME
fi

# write one line at the end of OUT4_FNAME since the next command skips last line
echo " " >> $OUT4_FNAME
echo " " >> $OUT4_FNAME

OUT_TAG=5
python3 _extract_pattern.py $OUT4_FNAME "(?<=\n).+?(?=\n)" -outtag="$OUT_TAG"\
    -offset=5 -prefix="ConvOfflineCheetah(conv_ntt, " -suffix=");"  --flatten

OUT5_FNAME=$(generate_out_fname "$OUT4_FNAME" "$OUT_TAG")

if [ $KEEP_INTERMEDIATE_FILES -eq 0 ]; then
    rm $OUT4_FNAME
fi

mv $OUT5_FNAME $(generate_out_fname "$INP_FNAME" "conv_off")