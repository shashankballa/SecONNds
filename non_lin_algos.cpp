
void funcMaxPool(int outSize, // output size
                int wSize, // flattened window size
                int bw, // data bitwidth
                int *inpSecShare, // input secret share morphed to wSize x outSize such that each column is a window
                int *outSecShare // flattenned output secret share
                ) {
    int *max_temp = new int[outSize];
    for (int r = 0; r < outSize; r++) {
        max_temp[r] = inpSecShare[r * wSize];
    }
    for (int c = 1; c < wSize; c++) {
        int *compare_with = new int[outSize];
        for (int r = 0; r < outSize; r++) {
            compare_with[r] = max_temp[r] - inpSecShare[r * wSize + c];
        }
        relu_oracle->relu(max_temp, compare_with, outSize);
        for (int r = 0; r < outSize; r++) {
            max_temp[r] += inpSecShare[r * wSize + c];
        }
    }
    for (int r = 0; r < outSize; r++) {
        outSecShare[r] = max_temp[r];
        outSecShare[r] &= ((1ULL << bw) - 1);
    }
}