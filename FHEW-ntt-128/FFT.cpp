#include <complex.h>
#include <fftw3.h>
#include "FFT.h"
#include "params.h"
#include "LWE.h"
#include <iostream>


double *in;
fftw_complex *out;
fftw_plan plan_fft_forw, plan_fft_back;
  
void FFTsetup() {
  in = (double*) fftw_malloc(sizeof(double) * 2*N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N + 2));
  plan_fft_forw = fftw_plan_dft_r2c_1d(2*N, in, out,  FFTW_PATIENT);
  plan_fft_back = fftw_plan_dft_c2r_1d(2*N, out, in,  FFTW_PATIENT);
}
  
void FFTforward(Ring_FFT res, const Ring_ModQ val) {
  for (int k = 0; k < N; ++k)	{
    in[k] = (double) (val[k]);
    in[k+N] = 0.0;			
  }
  fftw_execute(plan_fft_forw); 
  for (int k = 0; k < N2; ++k) 
    res[k] = (double complex) out[2*k+1];				
}

void FFTbackward(Ring_ModQ res, const Ring_FFT val){
  for (int k = 0; k < N2; ++k) {
    out[2*k+1] = (double complex) val[k]/N;
    out[2*k]   = (double complex) 0;
  }
  fftw_execute(plan_fft_back); 
  for (int k = 0; k < N; ++k)	
    res[k] = LWE::ModAddQ((ZmodQ) 0, (int)round(in[k]));
}

void AntiCyclicConv(Ring_ModQ val, const Ring_ModQ op1, const Ring_ModQ op2){
  //init phase
  Ring_ModQ tmp;
  ZmodQ extra;
  tmp[0] = op1[0];
  for (int i = 0; i < N; i++)
      val[i] = 0;
  for (int i = 1; i < N; i++) //[a[0], -a[N-1], -a[N-2], ... , -a[1]]
      tmp[N-i] = LWE::ModMulQ(op1[i], -1); 
  for (int i = 0; i < N; i++)
     //val[0] += tmp[i]*op2[i];
      val[0] = LWE::ModAddQ(val[0], LWE::ModMulQ(tmp[i], op2[i])); 

  //rotation phase
  for (int i = 1; i < N; i++){
    extra = tmp[N-1];
    for (int j = N-1; j > 0; j--) //anti-right rotate
      tmp[j] = tmp[j-1];
    tmp[0] = LWE::ModMulQ(extra, -1);

    for (int j = 0; j < N; j++)
      //val[i] += tmp[j]*op2[j]; 
      val[i] = LWE::ModAddQ(val[i], LWE::ModMulQ(tmp[j], op2[j]));     
  }
}


// a simple NTT implementation follows Nayuki-FFT-Portable

// Private data structure 
typedef struct 
{
  size_t n;
  size_t *bit_reversed;
  ZmodQ *omega_table;
  ZmodQ *omega_inv_table;
  ZmodQ *psi_table;
  ZmodQ *psi_inv_table;
} NttTables;

//declare structure
NttTables *tbl;

// Private function prototypes
static int32_t floor_log2(size_t n);
static size_t reverse_bits(size_t x, size_t n);

//Return a pointer to a structure of NTT tables. n must be a power of 2.
void NTTsetup( ){
  //Allocate tables
  tbl = (NttTables*) malloc(sizeof(NttTables));
  tbl->n = N;
  //Allocate arrays
  tbl->bit_reversed = (size_t*) malloc(N*sizeof(size_t));
  tbl->omega_table = (ZmodQ*) malloc(N*sizeof(ZmodQ));
  tbl->omega_inv_table = (ZmodQ*) malloc(N*sizeof(ZmodQ));
  tbl->psi_table = (ZmodQ*) malloc(N*sizeof(ZmodQ));
  tbl->psi_inv_table = (ZmodQ*) malloc(N*sizeof(ZmodQ));
  tbl->omega_table[0] = 1;
  tbl->omega_inv_table[0] = 1;
  tbl->psi_table[0] = 1;
  tbl->psi_inv_table[0] = 1;

  //Precompute values and store to tables
  int32_t levels = floor_log2(N);
  for (size_t i = 0; i < N; ++i)
    tbl->bit_reversed[i] = reverse_bits(i, levels);
  for (size_t i = 1; i < N; ++i){
    tbl->omega_table[i] = LWE::ModMulQ(tbl->omega_table[i-1], Omega);
    tbl->omega_inv_table[i] = LWE::ModMulQ(tbl->omega_inv_table[i-1], Omega_inverse);
    tbl->psi_table[i] = LWE::ModMulQ(tbl->psi_table[i-1], Psi);
    tbl->psi_inv_table[i] = LWE::ModMulQ(tbl->psi_inv_table[i-1], Psi_inverse);
  }
  // for (int i = 0; i < N; i++)
  //   std::cout << "bit_reversed[" << i << "] = " << tbl->bit_reversed[i] << "\n";
}

void NTTforward(Ring_NTT res, const Ring_ModQ val){
  size_t *bit_reversed = tbl->bit_reversed;
  size_t n = tbl->n;
  //Pre-process N points of val via Psi list
  ZmodQ *psi = tbl->psi_table;
  for (size_t i = 0; i < n; i++){
    res[i] = LWE::ModMulQ(val[i], psi[i]);
  }  
  //standard N-point-wise NTT forward transform
  // Bit-reversed addressing permutation
  for (size_t i = 0; i < n; i++){
    size_t j = bit_reversed[i];
    if (i < j){
      ZmodQ tp0 = res[i];
      ZmodQ tp1 = res[j];
      res[i] = tp1;
      res[j] = tp0; 
    }   
  }
  // Cooley-Tukey decimation-in-time radix-2 FFT
  ZmodQ *omega = tbl->omega_table;
  for (size_t size = 2; size <= n; size *= 2){
    size_t halfsize = size/2;
    size_t tablestep = n/size;
    for (size_t i = 0; i < n; i += size){
      size_t j, k;
      for (j = i, k = 0; j < i + halfsize; j++, k += tablestep){
        ZmodQ tp = LWE::ModMulQ(res[j+halfsize], omega[k]);
        res[j+halfsize] = LWE::ModSubQ(res[j], tp);
        res[j] = LWE::ModAddQ(res[j], tp);
      }
    }
    if (size == n) // Prevent overflow in 'size *= 2'
        break;
  }
}

void NTTbackward(Ring_ModQ res, const Ring_NTT val){
  size_t *bit_reversed = tbl->bit_reversed;
  size_t n = tbl->n;
  //standard N-point-wise NTT backward transform
  // Bit-reversed addressing permutation
  for (size_t i = 0; i < n; i++){
    size_t j = bit_reversed[i];
    res[i] = val[j];
  }
  // Cooley-Tukey decimation-in-time radix-2 FFT
  ZmodQ *omega = tbl->omega_inv_table;
  for (size_t size = 2; size <= n; size *= 2){
    size_t halfsize = size/2;
    size_t tablestep = n/size;
    for (size_t i = 0; i < n; i += size){
      size_t j, k;
      for (j = i, k = 0; j < i + halfsize; j++, k += tablestep){
        ZmodQ tp = LWE::ModMulQ(res[j+halfsize], omega[k]);
        res[j+halfsize] = LWE::ModSubQ(res[j], tp);
        res[j] = LWE::ModAddQ(res[j], tp);
      }
    }
    if (size == n) // Prevent overflow in 'size *= 2'
        break;
  }
  for (size_t i = 0; i < n; i++)
    res[i] = LWE::ModMulQ(res[i], N_inv); 
  //Post-process N points of val via Psi_inv list
  ZmodQ *psi = tbl->psi_inv_table;
  for (size_t i = 0; i < n; i++){
    res[i] = LWE::ModMulQ(res[i], psi[i]);
  }  
}

//Deallocate the given structure of  NTT tables
void NTTdestroy(){
  if (tbl == NULL)
    return;
  free(tbl->bit_reversed);
  free(tbl->omega_table);
  free(tbl->omega_inv_table);
  free(tbl->psi_table);
  free(tbl->psi_inv_table);
  free(tbl);
}

static int32_t floor_log2(size_t n){ // return the largest i such that 2^i <= n
  int32_t result = 0;
  for(; n > 1; n /= 2)
    result++;
  return result;
}

static size_t reverse_bits(size_t x, size_t n){ //return the bit reversal of the n-bit unsigned integer x
  size_t result = 0;
  for (size_t i = 0; i < n; i++, x >>= 1)
    result = (result << 1) | (x & 1);
  return result;
} 
