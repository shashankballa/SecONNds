
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
    -prefix="vector<vector<vector<vector<vector<seal::Plaintext>>>>> " -suffix="_pts;" -offset=6
OUT1_FNAME=$(generate_out_fname "$OUT0_FNAME" "$OUT_TAG")

mv $OUT1_FNAME $(generate_out_fname "$INP_FNAME" "tmp_pts")

OUT_TAG=2
python3 _extract_pattern.py $OUT0_FNAME "Conv2DWrapper\(([^,]+(?:, [^,]+)*?), [^,]+\);" -outtag="$OUT_TAG"\
    -offset=6

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

OUT_TAG=5
python3 _modify_and_insert.py $OUT4_FNAME 13 14 -outtag="$OUT_TAG" \
    -offset=6 -prefix="secret_share_vec_" -suffix="" -dropf=3 
OUT5_FNAME=$(generate_out_fname "$OUT4_FNAME" "$OUT_TAG")

if [ $KEEP_INTERMEDIATE_FILES -eq 0 ]; then
    rm $OUT4_FNAME
fi

OUT_TAG=6
python3 _modify_and_insert.py $OUT5_FNAME 13 14 -outtag="$OUT_TAG" \
    -offset=6 -prefix="noise_pts_" -suffix="" -dropf=3 
OUT6_FNAME=$(generate_out_fname "$OUT5_FNAME" "$OUT_TAG")

if [ $KEEP_INTERMEDIATE_FILES -eq 0 ]; then
    rm $OUT5_FNAME
fi

OUT_TAG=7
python3 _modify_and_insert.py $OUT6_FNAME 13 14 -outtag="$OUT_TAG" \
    -offset=6 -prefix="noise_cts_" -suffix="" -dropf=3 
OUT7_FNAME=$(generate_out_fname "$OUT6_FNAME" "$OUT_TAG")

if [ $KEEP_INTERMEDIATE_FILES -eq 0 ]; then
    rm $OUT6_FNAME
fi

# write one line at the end of OUT4_FNAME since the next command skips last line
echo " " >> $OUT7_FNAME
echo " " >> $OUT7_FNAME

OUT_TAG=8
python3 _extract_pattern.py $OUT7_FNAME "(?<=\n).+?(?=\n)" -outtag="$OUT_TAG"\
    -offset=5 -prefix="ConvOfflineHeliks(conv_ntt, " -suffix=");"  --flatten

OUT8_FNAME=$(generate_out_fname "$OUT7_FNAME" "$OUT_TAG")

if [ $KEEP_INTERMEDIATE_FILES -eq 0 ]; then
    rm $OUT7_FNAME
fi

mv $OUT8_FNAME $(generate_out_fname "$INP_FNAME" "conv_off")