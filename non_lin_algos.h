void send(int *data){};
void send(int data){};

void receive(int *data){};
void receive(int data){};


// y = sel * x
void multiplexer(
    // selection bits
    int sel,
    // input vector
    int x,
    // output vector
    int y,
    // bitwidth of output
    int b){

  // delta for COT
  int d = (x * (1 - 2 * sel)) % (1 << b);

  int m_s; // sender output message for COT, m0
  int m_r; // receiver output message for COT, m1
  if (party == 0) {
    send_cot(m_s, d, b);
    recv_cot(m_r, sel, b);
  } else {  // party == 1
    recv_cot(m_r, sel, b);
    send_cot(m_s, d, b);
  }

  y = ((x * sel + m_r - m_s) % (1 << b));
};


// Simple MSB to Wrap computation
void msb1_to_wrap(
    // input
    int x,
    // output share of Wrap(x)
    int w,
    // bitwidth of x
    int b) {
    
  int msb_xb = (x >> (b - 1)) & 1;
      
  if (party == 0) {
    
    w = prg.get(bool);

    int *m_s = new int[2]; // sender's OT messages
    m_s[0] = w;
    m_s[1] = w ^ msb_xb;
    sendOT(m_s);

  } else {  // party == 1

    bool c = msb_xb; // receiver's OT choice bit
    recvOT(w, c, 1);

  }
}

// Boolean to Arithmetic Shares
void B2A(
    // input boolean secret share
    int x,
    // output arithmetic secret share
    int y,
    // bitwidth of y
    int b){

  int N = (1 << b);

  if (party == 0) {
    int d = (-2 * x) & N; // delta for COT
    int m_s = send_cot(d, b);
    y = (x - m_s) % N;

  } else {  // party == 1

    int c = x; // choice for COT
    int m_r = recvOT(m_r, c, 1);
    y = (x + m_r) % N;
  }
}

void sec_and(
  // input secret shares
  int i0SS, 
  int i1SS){

  {aSS, bSS, cSS} = tripleGen.get(1);

  // Compute shares of correction bits for Secure AND with triples
  int eSS = aSS ^ i0SS
  int fSS = bSS ^ i1SS;

  // Communicate the shares of correction bits
  int eSS_;
  int fSS_;
  if (party == 0) {
    send(eSS);
    send(fSS);
    receive(eSS_);
    receive(fSS_);
  } else { // party = 1
    receive(eSS_);
    receive(fSS_);
    send(eSS);
    send(fSS);
  }

  e = eSS ^ eSS_;
  f = fSS ^ fSS_;

  // Secret shares of less than results for each bit
  int bSS_cmp;

  // Conpute the output shares for secure AND with triples

  if (party == 0)
    bSS_cmp = e & f;
  else // party = 1
    bSS_cmp = 0;

  bSS_cmp = bSS_cmp ^ (e & bSS) ^ (f & aSS) ^ cSS; 
}

void millionaire_compare(
      int oSS, int inp, int b,
      int party = 0,
      bool greater_than = true){

    {aSS, bSS, cSS} = tripleGen.get(2*b-1)

    int *i0SS = new int [b]; // secret shares of Alice's input
    int *i1SS = new int [b]; // secret shares of Bob's input
    std::fill_n(i0SS, b, 0);
    std::fill_n(i1SS, b, 0);
    
    for (int i = 0; i < b; i++) {
      if (party == 0) {
          i0SS[i] = (inp >> i) & 1;
          if(!greater_than) i0SS[i] ^= 1;
      } else { // party = 1
          i1SS[i] = (inp >> i) & 1;
          if(greater_than) i1SS[i] ^= 1;
      }
    }

    int *bSS_eql = new int[b];
    for (int i = 0; i < b; i++) {
        bSS_eql[i] = i0SS[i] ^ i1SS[i];
    }

    int *eSS = new int[b];
    int *fSS = new int[b];

    // Compute shares of correction bits for Secure AND with triples
    for (int i = 0; i < b; i ++) {
      eSS[i] = aSS[i]^ i0SS[i];
      fSS[i] = bSS[i]^ i1SS[i];
    }

    // Communicate the shares of correction bits
    int *e  = new int[b];
    int *f  = new int[b];
    if (party == 0) {
      send(eSS);
      send(fSS);
      receive(e);
      receive(f);
    } else { // party = 1
      receive(e);
      receive(f);
      send(eSS);
      send(fSS);
    }

    // Reconstruct the correction bits
    for (int i = 0; i < b; i++) {
      e[i] ^= eSS[i];
      f[i] ^= fSS[i];
    }

    // Secret shares of less than results for each bit
    int *bSS_cmp = new int[b];

    // Conpute the output shares for secure AND with triples
    for (int i = 0; i < b; i++) {

      if (party == 0)
        bSS_cmp[i] = e[i] & f[i];
      else // party = 1
        bSS_cmp[i] = 0;

      bSS_cmp[i] = bSS_cmp[i] ^ (e[i] & bSS[i]) ^ (f[i] & aSS[i]) ^ cSS[i];
    }
    
    for(int i = 0; i < b - 1; i++){
      
      // Compute the correction bits for Secure AND with triples
      int eSS_ = aSS[b + i] ^ bSS_cmp[i];
      int fSS_ = bSS[b + i] ^ bSS_eql[i+1];

      // Communicate the correction bits
      int e_, f_; 
      if (party == 0) {
        send(eSS_);
        send(fSS_);
        receive(e_ );
        receive(f_ );
      } else // party = 1
      {
        receive(e_ );
        receive(f_ );
        send(eSS_);
        send(fSS_);
      }       
      
      e_ ^= eSS_;
      f_ ^= fSS_;

      int eq_and_lt;
      if (party == 0)
        eq_and_lt = e_ & f_;
      else // party = 1
        eq_and_lt = 0;
      eq_and_lt = eq_and_lt ^ (e_ & bSS[b + i]) ^ (f_ & aSS[b + i]) ^ cSS[b + i];

      bSS_cmp[i + 1] ^= eq_and_lt;
    }

    oSS = bSS_cmp[b - 1];
};

void relu(int outSS, int inpSS, int b = 32, int party = 0,
        bool do_trunc = false) {

    int max_pos = (1 << (b - 1)) - 1; // max +ve number in ring mod 2^b
    int inpSS_msb = (inpSS >> (b - 1));
    int inpSS_abs = inpSS & max_pos;

    int inpSS_mill;
    if(party == 0) {
      inpSS_mill = inpSS_abs;
    } else { // party == 1
      inpSS_mill = max_pos - inpSS_abs; // This value is never negative.
    }

    // DRelu(x) = 1 ^ msb(x)
    //          = 1 ^ msb(x0) ^ msb(x1) ^ wrap
    // wrap     = 1(x0 + x1 > {2^{b-1} - 1})
    int wrapSS;
    millionaire_compare(wrapSS, inpSS_mill, b - 1);

    inpSS_msb = (inpSS_msb + wrapSS) % 2;
    if (party == 0) {
      inpSS_msb = inpSS_msb ^ 1;
    }
    multiplexer(inpSS_msb, inpSS, outSS, b, b);
};

void relu(int *outSS, int *inpSS, int b = 32, int party = 0,
        bool do_trunc = false) {};

void MaxPool(int outSize, // output size
                int wSize, // flattened window size
                int b, // data bitwidth
                int *inpSS, // input secret share morphed to wSize x outSize such that each column is a window
                int *outSS // flattenned output secret share
                ) {
  for (int r = 0; r < outSize; r++) {
    outSS[r] = inpSS[r * wSize];
  }
  for (int c = 1; c < wSize; c++) {
    int *diffSS = new int[outSize];
    for (int r = 0; r < outSize; r++) {
      diffSS[r] = outSS[r] - inpSS[r * wSize + c];
    }
    relu(outSS, diffSS, outSize);
    for (int r = 0; r < outSize; r++) {
      outSS[r] += inpSS[r * wSize + c];
    }
  }
};

// Truncate (right-shift) by shift in the same ring (round towards -inf)
void truncate(
      // Party ID
      int party,
      // input vector
      int iSS,
      // output vector
      int oSS,
      // right shift amount
      int s,
      // Input and output bitwidth
      int b){

  int N = (1 << b); // ring modulus

  // Ref "Secure evaluation of quantized neural networks"
  // https://eprint.iacr.org/2019/131.pdf
  if (party == 1) {
    const int m = b - 3;

    int big_positive = 1 << m;
    
    int adjust = (iSS + big_positive) % N;

    truncate_msb0(party, adjust, oSS, s, b);
    // Tr(x + 2^m, 2^f) = x/2^f + 2^{m - f}
    
    int offset = 1 << (m - s);
    oSS = (oSS - offset) % N;

  } else {
    truncate_msb0(party, iSS, oSS, s, b);
  }
}


// Truncate (right-shift) by shift in the same ring (round towards -inf)
// All elements have msb equal to 0.
void truncate_msb0(
    // party ID
    int party,
    // input secret share
    int iSS,
    // output vector
    int oSS,
    // right shift amount
    int s,
    // Input and output bitwidth
    int b) {

  int N = (1 << b); // ring modulus
  int N_s = N >> s; // truncated modulus

  int inSS_w;
  if (party == 0) {
    inSS_w = (iSS + N/2) % N;
  } else {
    inSS_w = iSS;
  }

  int inSS_msb = (inSS_w >> (b - 1)) & 1;
      
  int wSS_b;
  if (party == 0) {
    wSS_b = prg.get(bool);
    int *m_s = new int[2]; // sender's OT messages
    m_s[0] = wSS_b;
    m_s[1] = wSS_b ^ inSS_msb;
    sendOT(m_s);

  } else {  // party == 1

    bool c = inSS_msb; // receiver's OT choice bit
    wSS_b = recvOT(c, 1);
  }

  int wSS;
  if (party == 0) {
    int d = (-2 * wSS_b) & N; // delta for COT
    int m_s = sendCOT(d, b);
    wSS = (wSS_b - m_s) % N;

  } else {  // party == 1

    int c = wSS_b; // choice for COT
    int m_r = recvCOT(c, 1);
    
    wSS = (wSS_b + m_r) % N;
  }

  oSS = (((iSS >> s) % N_s) - N_s * wSS) % N;

  if (party == 0) {
      oSS = (oSS - N_s/2) % N;
  }
}


