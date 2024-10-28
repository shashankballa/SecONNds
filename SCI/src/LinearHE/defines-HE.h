/*
Authors: Deevashwer Rathee
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

#ifndef DEFINES_HE_H__
#define DEFINES_HE_H__

#include <cstdint>
#include <cmath>
#include <vector>

// #define HE_DEBUG 1

extern uint64_t prime_mod;
extern int32_t bitlength;
extern int32_t num_threads;

const uint64_t POLY_MOD_DEGREE_SMALL = 8192;
const uint64_t POLY_MOD_DEGREE = 8192;
const uint64_t POLY_MOD_DEGREE_LARGE = 32768;
const int32_t SMUDGING_BITLEN = 100 - bitlength;

// Default coefficient modulus for CryptFlow2
inline std::vector<int> GET_COEFF_MOD_CF2() { return {60, 60, 60, 49}; }

// Default coefficient modulus for HELiKs
inline std::vector<int> GET_COEFF_MOD_HLK() {
    
    std::vector<int> coeff_mod = {60, 60, 60};

    switch (bitlength) {
    case 37:
        coeff_mod = {55, 40, 60};
        break;
    case 32:
        coeff_mod = {43, 37, 50};
        break;
    }

    return coeff_mod;
}

/* Helper function for rounding to the next power of 2
 * Credit:
 * https://stackoverflow.com/questions/466204/rounding-up-to-next-power-of-2 */
inline int next_pow2(int val) { return pow(2, ceil(log(val) / log(2))); }

#endif // DEFINES_HE_H__
