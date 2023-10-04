#ifndef PARAM_H
#define PARAM_H

#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdint.h>

const int n = 512;

const int N = 1024;
const int N2 = N/2+1;

const int K = 4;
const int K2 = 2*K;

const int64_t Q = 134215681; // 2^27 - 2047, has 2N primitive root of unity
const int q = 512;

typedef int32_t ZmodQ;
typedef uint32_t uZmodQ;

const ZmodQ Psi = 4073518;
const ZmodQ Omega = 61607251;
const ZmodQ Psi_inverse = 18038264;
const ZmodQ Omega_inverse = 26901077;

const ZmodQ N_inv = -131070;  // N*N_inv = 1 mod Q
const ZmodQ v = 16776960;       // v = round(Q/2t) = round(Q/8)
const ZmodQ v_inverse = -8; // 1/v mod Q

// const int n = 512;

// const int N = 1024;
// const int N2 = N/2+1;

// const int K = 3;
// const int K2 = 2*K;

// const int64_t Q = 4294957057; // 2^32 - 10239, has 2N primitive root of unity
// const int q = 512;
// const int q2 = 256;

// typedef int32_t ZmodQ;
// typedef uint32_t uZmodQ;

// const ZmodQ Omega = 1991688061;
// const ZmodQ Omega_inverse = 28452606;
// const ZmodQ Psi = -1663203887;
// const ZmodQ Psi_inverse = 498062085;

// const ZmodQ N_inv = -4194294;  // N*N_inv = 1 mod Q
// const ZmodQ v = 536869632;       // v = round(Q/8)
// const ZmodQ v_inverse = -8; // 1/v mod Q

const ZmodQ vgprime[K] = {v, \
	(int64_t(v)<<7)%Q <= Q/2 ? (int64_t(v)<<7)%Q : (int64_t(v)<<7)%Q-Q, \
	(int64_t(v)<<14)%Q <= Q/2 ? (int64_t(v)<<14)%Q : (int64_t(v)<<14)%Q-Q, \
	(int64_t(v)<<21)%Q <= Q/2 ? (int64_t(v)<<21)%Q : (int64_t(v)<<21)%Q-Q};
const int g_bits[K] = {7, 7, 7, 7};
const int g_bits_31[K] = {25, 25, 25, 25}; // set g_bits + g_bits_31 = 32

/*const ZmodQ vgprime[K] = {v, \
	(int64_t(v)<<6)%Q <= Q/2 ? (int64_t(v)<<6)%Q : (int64_t(v)<<6)%Q-Q, \
	(int64_t(v)<<12)%Q <= Q/2 ? (int64_t(v)<<12)%Q : (int64_t(v)<<12)%Q-Q, \
	(int64_t(v)<<18)%Q <= Q/2 ? (int64_t(v)<<18)%Q : (int64_t(v)<<18)%Q-Q, \
	(int64_t(v)<<24)%Q <= Q/2 ? (int64_t(v)<<24)%Q : (int64_t(v)<<24)%Q-Q};
const int g_bits[K] = {6, 6, 6, 6, 6};
const int g_bits_31[K] = {26, 26, 26, 26, 26};*/


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
