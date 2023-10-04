#ifndef PARAM_H
#define PARAM_H

#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdint.h>

// const int n = 512;

// const int N = 1024;
// const int N2 = N/2+1;

// const int K = 4;
// const int K2 = 2*K;

// const int64_t Q = 134215681; // 2^27 - 2047, has 2N primitive root of unity
// const int q = 512;
// const int q2 = 256;

// typedef int32_t ZmodQ;
// typedef uint32_t uZmodQ;

// const ZmodQ Psi = 4073518;
// const ZmodQ Omega = 61607251;
// const ZmodQ Psi_inverse = 18038264;
// const ZmodQ Omega_inverse = 26901077;

// const ZmodQ N_inv = -131070;  // N*N_inv = 1 mod Q
// const ZmodQ v = 16776960;       // v = round(Q/2t) = round(Q/8)
// const ZmodQ v_inverse = -8; // 1/v mod Q

const int n = 1024;

const int N = 2048;
const int N2 = N/2+1;

const int K = 5;
const int K2 = 2*K;

const int64_t Q = 536813569; // 2^29-57343, has 2N primitive root of unity
const int q = 512;
const int q2 = q/2;

typedef int32_t ZmodQ;
typedef uint32_t uZmodQ;

const ZmodQ Psi = 251063407;
const ZmodQ Omega = -138213656;
const ZmodQ Psi_inverse = -108056907;
const ZmodQ Omega_inverse = -132913493;

const ZmodQ N_inv = -262116;  // N*N_inv = 1 mod Q
const ZmodQ v = 67101696;       // v = round(Q/2t) = round(Q/8)
const ZmodQ v_inverse = -8; // 1/v mod Q

/*const ZmodQ vgprime[K] = {v, \
	(int64_t(v)<<10)%Q <= Q/2 ? (int64_t(v)<<10)%Q : (int64_t(v)<<10)%Q-Q, \
	(int64_t(v)<<20)%Q <= Q/2 ? (int64_t(v)<<20)%Q : (int64_t(v)<<20)%Q-Q};
const int g_bits[K] = {10, 10, 10};
const int g_bits_31[K] = {21, 21, 21}; // set g_bits + g_bits_31 = 31*/


/*const ZmodQ vgprime[K] = {v, \
	(int64_t(v)<<8)%Q <= Q/2 ? (int64_t(v)<<8)%Q : (int64_t(v)<<8)%Q-Q, \
	(int64_t(v)<<16)%Q <= Q/2 ? (int64_t(v)<<16)%Q : (int64_t(v)<<16)%Q-Q, \
	(int64_t(v)<<24)%Q <= Q/2 ? (int64_t(v)<<24)%Q : (int64_t(v)<<24)%Q-Q};
const int g_bits[K] = {8, 8, 8, 8};
const int g_bits_31[K] = {24, 24, 24, 24};*/

const ZmodQ vgprime[K] = {v, \
	(int64_t(v)<<6)%Q <= Q/2 ? (int64_t(v)<<6)%Q : (int64_t(v)<<6)%Q-Q, \
	(int64_t(v)<<12)%Q <= Q/2 ? (int64_t(v)<<12)%Q : (int64_t(v)<<12)%Q-Q, \
	(int64_t(v)<<18)%Q <= Q/2 ? (int64_t(v)<<18)%Q : (int64_t(v)<<18)%Q-Q, \
	(int64_t(v)<<24)%Q <= Q/2 ? (int64_t(v)<<24)%Q : (int64_t(v)<<24)%Q-Q};
const int g_bits[K] = {6, 6, 6, 6, 6};
const int g_bits_31[K] = {26, 26, 26, 26, 26};

const int KS_base = 25;
const int KS_exp = 7;
const ZmodQ KS_table[KS_exp] = {1,
			   25,
			   25*25,
			   25*25*25,
			   25*25*25*25,
			   25*25*25*25*25,
			   25*25*25*25*25*25};

const int BS_base = 23;
const int BS_exp = 2;
const int BS_table[BS_exp] = {1,23};
// const int BS_base = 2;
// const int BS_exp = 9;
// const int BS_table[9] = {1,2,4,8,16,32,64,128,256};

typedef ZmodQ Ring_ModQ[N];
typedef fftw_complex Ring_FFT[N2];
typedef ZmodQ Ring_NTT[N];

enum BinGate {OR, AND, NOR, NAND};
const int GateConst[4] = {15*q/8, 9*q/8, 11*q/8, 13*q/8};


#endif
