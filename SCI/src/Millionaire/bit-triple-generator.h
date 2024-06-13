/*
Authors: Deevashwer Rathee
Modified by Zhicong Huang.
Copyright:
Copyright (c) 2020 Microsoft Research
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef TRIPLE_GENERATOR_H__
#define TRIPLE_GENERATOR_H__
#include "OT/emp-ot.h"

#define PRINT_TIME 1
#define PRINT_COMP 0
#define PRINT_COMM 1
#if PRINT_COMP || PRINT_COMM || PRINT_TIME
#include <iomanip>
#endif

enum TripleGenMethod {
  Ideal,          // (Insecure) Ideal Functionality
  _2ROT,          // 1 Bit Triple from 2 ROT
  _16KKOT_to_4OT, // 2 Bit Triples from 1oo16 KKOT to 1oo4 OT
  _8KKOT,         // 2 Correlated Bit Triples from 1oo8 KKOT
};

class Triple {
public:
  bool packed;
  uint8_t *ai;
  uint8_t *bi;
  uint8_t *ci;
  int num_triples, num_bytes, offset;

  Triple(int num_triples, bool packed = false, int offset = 0) {
    assert((offset < num_triples) || (num_triples == 0));
    this->num_triples = num_triples;
    this->packed = packed;
    if (packed) {
      assert(num_triples % 8 == 0);
      assert(offset % 8 == 0);
      this->num_bytes = num_triples / 8;
    } else
      this->num_bytes = num_triples;
    if (offset == 0)
      this->offset = 1;
    else
      this->offset = offset;
    assert((num_triples % this->offset) == 0);
    this->ai = new uint8_t[num_bytes];
    this->bi = new uint8_t[num_bytes];
    this->ci = new uint8_t[num_bytes];
  }

  ~Triple() {
    delete[] ai;
    delete[] bi;
    delete[] ci;
  }
};

#define BSIZE 67108864 // Default buffer size
#define CSIZE 67108864 // Default chunk size
template <typename IO> class TripleGenerator {
public:
  IO *io = nullptr;
  sci::OTPack<IO> *otpack = nullptr;
  sci::PRG128 *prg;
  int party;

  // Buffer implementation to use pre-generated triples from the offline phase
  uint8_t *Bai;         // Buffer for Ai
  uint8_t *Bbi;         // Buffer for Bi
  uint8_t *Bci;         // Buffer for Ci
  int buffSize  = (BSIZE >> 3) << 3; // Number of triples in the buffer (always multiple of 8)
  int buffBytes = buffSize >> 3;     // Number of bytes in the buffer (always = buffSize/8)
  int chunkSize = (CSIZE >> 3) << 3; // Number of triples to be generated in one go (always multiple of 8)
  bool buffEnable = false; // Flag to enable/disable buffer
  int buffPtr = 0;     // Pointer to the current triple in the buffer

  TripleGenerator(int party, IO *io, sci::OTPack<IO> *otpack
                  , bool enable_buffer = false
                  ) {
    this->io = io;
    this->otpack = otpack;
    this->prg = new sci::PRG128;

    if(enable_buffer) {

#if PRINT_COMM || PRINT_TIME
    std::string f_tag = "NEW | 3Gen";
#endif
#if PRINT_COMM
    int _wcomm = 10;
    std::stringstream log_comm;
    uint64_t start_comm = io->counter;
    uint64_t total_comm = 0;
    log_comm << "P" << party << " COMM | ";
    log_comm << f_tag;
    log_comm << ": enable_buffer = " << std::boolalpha << enable_buffer;
    log_comm << ", buffSize = " << buffSize;
    log_comm << ", chunkSize = " << chunkSize;
    log_comm << std::endl;
#endif
#if PRINT_TIME
    int _wtime = 10;
    std::stringstream log_time;
    auto start = std::chrono::system_clock::now();
    auto end   = std::chrono::system_clock::now();
    std::chrono::duration<double> total_time = end - start;
    log_time << "P" << party << " TIME | ";
    log_time << f_tag;
    log_time << ": enable_buffer = " << std::boolalpha << enable_buffer;
    log_time << ", buffSize = " << buffSize;
    log_time << ", chunkSize = " << chunkSize;
    log_time << std::endl;
#endif
      refillBuffer(party, buffSize);

#if PRINT_TIME
    // get running time of leaf OTs in ms
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> refill_time = end - (start + total_time);
    total_time += refill_time;
    log_time << "P" << party << " TIME | ";
    log_time << f_tag;
    log_time << ": refill = " << std::setw(_wtime) << refill_time.count() * 1000;
    log_time << " ms";
    log_time << std::endl;
#endif
#if PRINT_COMM
    uint64_t refill_comm = io->counter - (start_comm + total_comm);
    total_comm += refill_comm;
    log_comm << "P" << party << " COMM | ";
    log_comm << f_tag;
    log_comm << ": refill = " << std::setw(_wcomm) << refill_comm;
    log_comm << " bytes";
    log_comm << std::endl;
#endif
#if PRINT_COMM
    std::cout << log_comm.str();
#endif
#if PRINT_TIME
    std::cout << log_time.str();
#endif
    }
  }

  ~TripleGenerator() { 
    delete prg; 
    // delete buffer
    if(buffEnable) {
      delete[] Bai;
      delete[] Bbi;
      delete[] Bci;
    }
  }

  void generate(int party, uint8_t *ai, uint8_t *bi, uint8_t *ci,
                int num_triples, TripleGenMethod method, bool packed = false,
                int offset = 1) {
    if (!num_triples)
      return;
    switch (method) {
    case Ideal: {
      int num_bytes = ceil((double)num_triples / 8);
      if (party == sci::ALICE) {
        uint8_t *a = new uint8_t[num_bytes];
        uint8_t *b = new uint8_t[num_bytes];
        uint8_t *c = new uint8_t[num_bytes];
        if (packed) {
          prg->random_data(ai, num_bytes);
          prg->random_data(bi, num_bytes);
          prg->random_data(ci, num_bytes);
        } else {
          prg->random_bool((bool *)ai, num_triples);
          prg->random_bool((bool *)bi, num_triples);
          prg->random_bool((bool *)ci, num_triples);
        }
        prg->random_data(a, num_bytes);
        prg->random_data(b, num_bytes);
        for (int i = 0; i < num_triples; i += 8) {
          c[i / 8] = a[i / 8] & b[i / 8];
          if (packed) {
            a[i / 8] ^= ai[i / 8];
            b[i / 8] ^= bi[i / 8];
            c[i / 8] ^= ci[i / 8];
          } else {
            uint8_t temp_a, temp_b, temp_c;
            if (num_triples - i >= 8) {
              temp_a = sci::bool_to_uint8(ai + i, 8);
              temp_b = sci::bool_to_uint8(bi + i, 8);
              temp_c = sci::bool_to_uint8(ci + i, 8);
            } else {
              temp_a = sci::bool_to_uint8(ai + i, num_triples - i);
              temp_b = sci::bool_to_uint8(bi + i, num_triples - i);
              temp_c = sci::bool_to_uint8(ci + i, num_triples - i);
            }
            a[i / 8] ^= temp_a;
            b[i / 8] ^= temp_b;
            c[i / 8] ^= temp_c;
          }
        }
        io->send_data(a, num_bytes);
        io->send_data(b, num_bytes);
        io->send_data(c, num_bytes);
        delete[] a;
        delete[] b;
        delete[] c;
      } else {
        if (packed) {
          io->recv_data(ai, num_bytes);
          io->recv_data(bi, num_bytes);
          io->recv_data(ci, num_bytes);
        } else {
          uint8_t *a = new uint8_t[num_bytes];
          uint8_t *b = new uint8_t[num_bytes];
          uint8_t *c = new uint8_t[num_bytes];
          io->recv_data(a, num_bytes);
          io->recv_data(b, num_bytes);
          io->recv_data(c, num_bytes);

          for (int i = 0; i < num_triples; i += 8) {
            if (num_triples - i >= 8) {
              sci::uint8_to_bool(ai + i, a[i / 8], 8);
              sci::uint8_to_bool(bi + i, b[i / 8], 8);
              sci::uint8_to_bool(ci + i, c[i / 8], 8);
            } else {
              sci::uint8_to_bool(ai + i, a[i / 8], num_triples - i);
              sci::uint8_to_bool(bi + i, b[i / 8], num_triples - i);
              sci::uint8_to_bool(ci + i, c[i / 8], num_triples - i);
            }
          }
          delete[] a;
          delete[] b;
          delete[] c;
        }
      }
      break;
    }
    case _2ROT: {
#if USE_CHEETAH
        uint8_t *a, *b, *c;
        if (packed) {
            a = new uint8_t[num_triples];
            b = new uint8_t[num_triples];
            c = new uint8_t[num_triples];
        } else {
            a = ai;
            b = bi;
            c = ci;
        }
        uint8_t *u, *v;
        u = new uint8_t[num_triples];
        v = new uint8_t[num_triples];
        io->sync();
        switch (party) {
            case sci::ALICE: {
                otpack->silent_ot_reversed->template recv_ot_rm_rc<uint8_t>(u, (bool*)a, num_triples, 1);
                otpack->silent_ot->send_ot_rm_rc(v, b, num_triples, 1);
                break;
            }
            case sci::BOB: {
                otpack->silent_ot_reversed->template send_ot_rm_rc<uint8_t>(v, b, num_triples, 1);
                otpack->silent_ot->recv_ot_rm_rc(u, (bool*)a, num_triples, 1);
                break;
            }
        }
        io->flush();

        for (int i = 0; i < num_triples; i++)
            b[i] = b[i] ^ v[i];
        for (int i = 0; i < num_triples; i++)
            c[i] = (a[i] & b[i]) ^ u[i] ^ v[i];


        delete[] u;
        delete[] v;
        if (packed) {
            for (int i = 0; i < num_triples; i += 8) {
                ai[i / 8] = sci::bool_to_uint8(a + i, 8);
                bi[i / 8] = sci::bool_to_uint8(b + i, 8);
                ci[i / 8] = sci::bool_to_uint8(c + i, 8);
            }
            delete[] a;
            delete[] b;
            delete[] c;
        }
#else
      throw std::invalid_argument("To be implemented");
#endif
      break;
    }
    case _16KKOT_to_4OT: {
      assert((num_triples & 1) == 0); // num_triples is even
      uint8_t *a, *b, *c;
      if (packed) {
        a = new uint8_t[num_triples];
        b = new uint8_t[num_triples];
        c = new uint8_t[num_triples];
      } else {
        a = ai;
        b = bi;
        c = ci;
      }
      prg->random_bool((bool *)a, num_triples);
      prg->random_bool((bool *)b, num_triples);
      switch (party) {
      case sci::ALICE: {
        prg->random_bool((bool *)c, num_triples);
        uint8_t **ot_messages; // (num_triples/2) X 16
        ot_messages = new uint8_t *[num_triples / 2];
        for (int i = 0; i < num_triples; i += 2)
          ot_messages[i / 2] = new uint8_t[16];
        for (int j = 0; j < 16; j++) {
          uint8_t bits_j[4]; // a01 || b01 || a11 || b11 (LSB->MSB)
          sci::uint8_to_bool(bits_j, j, 4);
          for (int i = 0; i < num_triples; i += 2) {
            ot_messages[i / 2][j] =
                ((((a[i + 1] ^ bits_j[2]) & (b[i + 1] ^ bits_j[3])) ^ c[i + 1])
                 << 1) |
                (((a[i] ^ bits_j[0]) & (b[i] ^ bits_j[1])) ^ c[i]);
          }
        }
        // otpack->kkot_16->send(ot_messages, num_triples/2, 2);
        otpack->kkot[3]->send(ot_messages, num_triples / 2, 2);
        for (int i = 0; i < num_triples; i += 2)
          delete[] ot_messages[i / 2];
        delete[] ot_messages;
        break;
      }
      case sci::BOB: {
        uint8_t *ot_selection = new uint8_t[(size_t)num_triples / 2];
        uint8_t *ot_result = new uint8_t[(size_t)num_triples / 2];
        for (int i = 0; i < num_triples; i += 2) {
          ot_selection[i / 2] =
              (b[i + 1] << 3) | (a[i + 1] << 2) | (b[i] << 1) | a[i];
        }
        // otpack->kkot_16->recv(ot_result, ot_selection, num_triples/2, 2);
        otpack->kkot[3]->recv(ot_result, ot_selection, num_triples / 2, 2);
        for (int i = 0; i < num_triples; i += 2) {
          c[i] = ot_result[i / 2] & 1;
          c[i + 1] = ot_result[i / 2] >> 1;
        }
        delete[] ot_selection;
        delete[] ot_result;
        break;
      }
      }
      if (packed) {
        for (int i = 0; i < num_triples; i += 8) {
          ai[i / 8] = sci::bool_to_uint8(a + i, 8);
          bi[i / 8] = sci::bool_to_uint8(b + i, 8);
          ci[i / 8] = sci::bool_to_uint8(c + i, 8);
        }
        delete[] a;
        delete[] b;
        delete[] c;
      }
      break;
    }
    case _8KKOT: {
      assert((num_triples & 1) == 0); // num_triples is even
      uint8_t *a, *b, *c;
      if (packed) {
        a = new uint8_t[num_triples];
        b = new uint8_t[num_triples];
        c = new uint8_t[num_triples];
      } else {
        a = ai;
        b = bi;
        c = ci;
      }
      for (int i = 0; i < num_triples; i += 2 * offset) {
        prg->random_bool((bool *)a + i, offset);
        memcpy(a + i + offset, a + i, offset);
      }
      prg->random_bool((bool *)b, num_triples);
      switch (party) {
      case sci::ALICE: {
        prg->random_bool((bool *)c, num_triples);
        uint8_t **ot_messages; // (num_triples/2) X 8
        ot_messages = new uint8_t *[num_triples / 2];
        for (int i = 0; i < num_triples; i += 2)
          ot_messages[i / 2] = new uint8_t[8];
        for (int j = 0; j < 8; j++) {
          uint8_t bits_j[3]; // a01 || b01 || b11 (LSB->MSB)
          sci::uint8_to_bool(bits_j, j, 3);
          for (int i = 0; i < num_triples; i += 2 * offset) {
            for (int k = 0; k < offset; k++) {
              ot_messages[i / 2 + k][j] =
                  ((((a[i + k] ^ bits_j[0]) & (b[i + offset + k] ^ bits_j[2])) ^
                    c[i + offset + k])
                   << 1) |
                  (((a[i + k] ^ bits_j[0]) & (b[i + k] ^ bits_j[1])) ^
                   c[i + k]);
            }
          }
        }
        // otpack->kkot_8->send(ot_messages, num_triples/2, 2);
        otpack->kkot[2]->send(ot_messages, num_triples / 2, 2);
        for (int i = 0; i < num_triples; i += 2)
          delete[] ot_messages[i / 2];
        delete[] ot_messages;
        break;
      }
      case sci::BOB: {
        uint8_t *ot_selection = new uint8_t[(size_t)num_triples / 2];
        uint8_t *ot_result = new uint8_t[(size_t)num_triples / 2];
        for (int i = 0; i < num_triples; i += 2 * offset) {
          for (int k = 0; k < offset; k++)
            ot_selection[i / 2 + k] =
                (b[i + offset + k] << 2) | (b[i + k] << 1) | a[i + k];
        }
        // otpack->kkot_8->recv(ot_result, ot_selection, num_triples/2, 2);
        otpack->kkot[2]->recv(ot_result, ot_selection, num_triples / 2, 2);
        for (int i = 0; i < num_triples; i += 2 * offset) {
          for (int k = 0; k < offset; k++) {
            c[i + k] = ot_result[i / 2 + k] & 1;
            c[i + offset + k] = ot_result[i / 2 + k] >> 1;
          }
        }
        delete[] ot_selection;
        delete[] ot_result;
        break;
      }
      }
      if (packed) {
        for (int i = 0; i < num_triples; i += 8) {
          ai[i / 8] = sci::bool_to_uint8(a + i, 8);
          bi[i / 8] = sci::bool_to_uint8(b + i, 8);
          ci[i / 8] = sci::bool_to_uint8(c + i, 8);
        }
        delete[] a;
        delete[] b;
        delete[] c;
      }
      break;
    }
    }
  }

  void generate(int party, Triple *triples, TripleGenMethod method) {
    generate(party, triples->ai, triples->bi, triples->ci, triples->num_triples,
             method, triples->packed, triples->offset);
  }

  /* Refill the buffer with new triples
    Resets the buffer size to num_triples
    Resets the chunk size to min(buffSize, CSIZE)
    Generates triples in chunks of size chunkSize with the specified method
    Resets the buffer pointer to 0
    Enables the buffer
  */
  void refillBuffer(int party, int num_triples, 
#if USE_CHEETAH
    TripleGenMethod method = _2ROT
#else
    TripleGenMethod method = _16KKOT_to_4OT
#endif
    ) {
    buffSize = num_triples;
    buffBytes = ceil((double)buffSize / 8);
    Bai = new uint8_t[buffBytes];
    Bbi = new uint8_t[buffBytes];
    Bci = new uint8_t[buffBytes];

    chunkSize = std::min(buffSize, chunkSize);
    int nChunks = ceil((double)buffSize / chunkSize);
    for (int i = 0; i < nChunks; i++) {
      int start = i * chunkSize;
      int startBytes = start / 8;
      int end = std::min(buffSize, (i + 1) * chunkSize);
      generate(party, 
                Bai + startBytes, 
                Bbi + startBytes, 
                Bci + startBytes,
                end - start, method, true);
    }
    buffPtr = 0;
    buffEnable  = true;
  }

  /* get triples from the buffer
    - If the buffer is disabled or triples are not packed, generates triples using the specified method
    - Fills the input triples object with triples from the buffer 
    - If the pointer reaches the end of the buffer, resets the buffer size to max(buffSize, num_triples) and refills
    - Copies triples from the buffer to the input triples object and increments the buffer pointer
  */
  void get(int party, Triple *triples, 
#if USE_CHEETAH
    TripleGenMethod method = _2ROT
#else
    TripleGenMethod method = _16KKOT_to_4OT
#endif
  ){
    if((!buffEnable) || (!triples->packed)){
      generate(party, triples->ai, triples->bi, triples->ci, triples->num_triples, method, triples->packed, triples->offset);
      return;
    }
    if(buffPtr+triples->num_triples > buffSize){
      int n_trips = triples->num_triples;
      int r_bytes = 0;
      if (buffPtr < buffSize){ // copy the remaining triples from the buffer
        r_bytes = (buffSize - buffPtr) >> 3;
        int startBytes = buffPtr >> 3;
        memcpy(triples->ai, Bai + startBytes, r_bytes);
        memcpy(triples->bi, Bbi + startBytes, r_bytes);
        memcpy(triples->ci, Bci + startBytes, r_bytes);
        buffPtr += (r_bytes << 3);
      }
      n_trips -= (r_bytes << 3);
      // generate triples in chunks of size chunkSize
      int nChunks = ceil((double)n_trips / chunkSize);
      for (int i = 0; i < nChunks; i++) {
        int start = i * chunkSize;
        int startBytes = start >> 3;
        int end = std::min(n_trips, (i + 1) * chunkSize);
        generate(party, 
                  triples->ai + r_bytes + startBytes, 
                  triples->bi + r_bytes + startBytes,
                  triples->ci + r_bytes + startBytes,
                  end - start, method, true);
      }
      if (buffPtr >= buffSize){
        refillBuffer(party, buffSize, method);
      }
      return;
    }
    if(buffPtr+triples->num_triples < buffSize){
      int startBytes = buffPtr >> 3;
      memcpy(triples->ai, Bai + startBytes, triples->num_bytes);
      memcpy(triples->bi, Bbi + startBytes, triples->num_bytes);
      memcpy(triples->ci, Bci + startBytes, triples->num_bytes);
      buffPtr += (triples->num_bytes) << 3;
      return;
    }
  }
};
#endif // TRIPLE_GENERATOR_H__
