#ifndef PARAM_H
#define PARAM_H

#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdint.h>

//#define _BARRET
const int n = 256;

const int N = 1024;
const int N2 = N/2+1;

const int K = 5;
const int K2 = 2*K;

const int64_t Q = 134215681; // 2^27 - 2047, has 2N primitive root of unity
const int q = 512;

const int barret_k = int( ceil(log(Q)/log(2)) ); 
const int64_t barret_r = ( int64_t(1)<<(2*barret_k) )/Q;



typedef int32_t ZmodQ;
typedef uint32_t uZmodQ;

typedef signed __int128 int128_t;
typedef unsigned __int128 uint128_t;

const ZmodQ Psi = 4073518;
const ZmodQ Omega = 61607251;
const ZmodQ Psi_inverse = 18038264;
const ZmodQ Omega_inverse = 26901077;

const ZmodQ N_inv = -131070;  // N*N_inv = 1 mod Q
const ZmodQ v = 16776960;       // v = round(Q/2t) = round(Q/8)
const ZmodQ v_inverse = -8; // 1/v mod Q


/*const ZmodQ vgprime[K] = {v, \
	(int64_t(v)<<7)%Q <= Q/2 ? (int64_t(v)<<7)%Q : (int64_t(v)<<7)%Q-Q, \
	(int64_t(v)<<14)%Q <= Q/2 ? (int64_t(v)<<14)%Q : (int64_t(v)<<14)%Q-Q, \
	(int64_t(v)<<21)%Q <= Q/2 ? (int64_t(v)<<21)%Q : (int64_t(v)<<21)%Q-Q};
const int g_bits[K] = {7, 7, 7, 7};
const int g_bits_31[K] = {25, 25, 25, 25}; // set g_bits + g_bits_31 = 32*/


const ZmodQ vgprime[K] = {v, \
	(int64_t(v)<<9)%Q <= Q/2 ? (int64_t(v)<<9)%Q : (int64_t(v)<<9)%Q-Q, \
	(int64_t(v)<<18)%Q <= Q/2 ? (int64_t(v)<<18)%Q : (int64_t(v)<<18)%Q-Q};
const int g_bits[K] = {9};
const int g_bits_31[K] = {23}; // set g_bits + g_bits_31 = 32, this parameter does not work


const int KS_base = 25;
const int KS_exp = 6;
const ZmodQ KS_table[KS_exp] = {1,
			   25,
			   25*25,
			   25*25*25,
			   25*25*25*25,
			   25*25*25*25*25};

const int BS_base = 23;
const int BS_exp = 2;
const int BS_table[2] = {1,23};
// const int BS_base = 2;
// const int BS_exp = 9;
// const int BS_table[9] = {1,2,4,8,16,32,64,128,256};

typedef ZmodQ Ring_ModQ[N];
typedef fftw_complex Ring_FFT[N2];
typedef ZmodQ Ring_NTT[N];

enum BinGate {OR, AND, NOR, NAND};
const int GateConst[4] = {15*q/8, 9*q/8, 11*q/8, 13*q/8};


#endif
