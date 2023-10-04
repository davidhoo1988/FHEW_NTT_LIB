#ifndef PARAM_H
#define PARAM_H

#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdint.h>


const int n = 512;

const int N = 2048;
const int N2 = N/2+1;

const int K = 3;
const int K2 = 2*K;

const int64_t Q = 137438822401; // 2^37-131071, has 2N primitive root of unity
const int q = 512;
//const int q2 = q/2;

typedef int64_t ZmodQ;
typedef uint64_t uZmodQ;

typedef signed __int128 int128_t;
typedef unsigned __int128 uint128_t;

const ZmodQ Psi = -39095648584;
const ZmodQ Omega = -61829332980;
const ZmodQ Psi_inverse = -20759711371;
const ZmodQ Omega_inverse = 13599156294;

const ZmodQ N_inv = -67108800;  // N*N_inv = 1 mod Q
const ZmodQ v = 17179852800;       // v = round(Q/2t) = round(Q/8)
const ZmodQ v_inverse = -8; // 1/v mod Q

const ZmodQ vgprime[K] = {v, \
	(int64_t(v)<<13)%Q <= Q/2 ? (int64_t(v)<<13)%Q : (int64_t(v)<<13)%Q-Q, \
	(int64_t(v)<<26)%Q <= Q/2 ? (int64_t(v)<<26)%Q : (int64_t(v)<<26)%Q-Q};
const int g_bits[K] = {13, 13, 13};
const int g_bits_31[K] = {51, 51, 51};

/*const ZmodQ vgprime[K] = {v, \
	(int64_t(v)<<10)%Q <= Q/2 ? (int64_t(v)<<10)%Q : (int64_t(v)<<10)%Q-Q, \
	(int64_t(v)<<20)%Q <= Q/2 ? (int64_t(v)<<20)%Q : (int64_t(v)<<20)%Q-Q, \
	(int128_t(v)<<30)%Q <= Q/2 ? (int128_t(v)<<30)%Q : (int128_t(v)<<30)%Q-Q};
const int g_bits[K] = {10, 10, 10, 10};
const int g_bits_31[K] = {53, 53, 53, 53};*/

/*const ZmodQ vgprime[K] = {v, \
	(int64_t(v)<<8)%Q <= Q/2 ? (int64_t(v)<<8)%Q : (int64_t(v)<<8)%Q-Q, \
	(int64_t(v)<<16)%Q <= Q/2 ? (int64_t(v)<<16)%Q : (int64_t(v)<<16)%Q-Q, \
	(int64_t(v)<<24)%Q <= Q/2 ? (int64_t(v)<<24)%Q : (int64_t(v)<<24)%Q-Q, \
	(int128_t(v)<<32)%Q <= Q/2 ? (int128_t(v)<<32)%Q : (int128_t(v)<<32)%Q-Q};
const int g_bits[K] = {8, 8, 8, 8, 8};
const int g_bits_31[K] = {55, 55, 55, 55, 55};*/


const int KS_base = 25;
const int KS_exp = 8; //KS_exp = ceil(log(Q, KS_base))
const ZmodQ KS_table[KS_exp] = {1,
			   25,
			   25*25,
			   25*25*25,
			   25*25*25*25,
			   25*25*25*25*25,
			   25*25*25*25*25*25,
			   int64_t(25)*25*25*25*25*25*25};

const int BS_base = 23;
const int BS_exp = 2; //BS_exp = ceil(log(q, BS_base))
const int BS_table[BS_exp] = {1, 23};
// const int BS_base = 2;
// const int BS_exp = 9;
// const int BS_table[9] = {1,2,4,8,16,32,64,128,256};

typedef ZmodQ Ring_ModQ[N];
typedef fftw_complex Ring_FFT[N2];
typedef ZmodQ Ring_NTT[N];

enum BinGate {OR, AND, NOR, NAND};
const int GateConst[4] = {15*q/8, 9*q/8, 11*q/8, 13*q/8};


#endif
