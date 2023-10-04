#include <iostream>
#include <cstdlib>
#include "FHEW.h"
#include "FFT.h"
#include <cassert>

using namespace std;
namespace FHEW {  
  Ring_FFT t_TestMSB;

  void Setup() {
    FFTsetup();
    Ring_ModQ tmsb;
    tmsb[0]=-1;
    for (int i = 1; i < N; ++i) tmsb[i]=1;
    FFTforward(t_TestMSB, tmsb);
  }
  
  void FHEWencrypt(ct_FFT ct, Ring_FFT sk_FFT, int m) { // Used by BootstrappingKeyGen
    Ring_FFT ai;
    ct_ModQ res;
    int mm = (((m % q) + q) % q) * (2*N/q);             // Reduce mod q (dealing with negative number as well)
    int sign = 1;
    if (mm >= N) { mm -= N; sign = -1; }
    
    for (int i = 0; i < K2; ++i) {
      for (int k = 0; k < N; ++k) res[i][0][k] = rand(); // % Q
      FFTforward(ai, res[i][0]);
      for (int k = 0; k < N2; ++k) 
        ai[k] = ((double complex) ai[k]) * ((double complex) sk_FFT[k]);
      FFTbackward(res[i][1], ai);
      for (int k = 0; k < N; ++k) res[i][1][k] += Sample(Chi1);    // Add error [a,as+e]
    }
    // Use Convolution directly
    // Ring_ModQ sk;
    // FFTbackward(sk, sk_FFT);
    // for (int i = 0; i < K2; ++i) {
    //   for (int k = 0; k < N; ++k) res[i][0][k] = rand(); // % Q
    //   AntiCyclicConv(res[i][1], res[i][0], sk);
    //   for (int k = 0; k < N; ++k) res[i][1][k] += Sample(Chi1);    // Add error [a,as+e]
    // }
    
    for (int i = 0; i < K; ++i) {
      res[2*i  ][0][mm] += sign*vgprime[i]; // Add G Multiple
      res[2*i+1][1][mm] += sign*vgprime[i]; // [a,as+e] + X^m *G
    }
    for (int i = 0; i < K2; ++i)
      for (int j = 0; j < 2; ++j)
        FFTforward(ct[i][j], res[i][j]);
  }
  
  int FHEWdecrypt(ct_FFT1 ct, Ring_FFT sk_FFT){
    Ring_FFT ai;
    Ring_ModQ asi;
    ct_ModQ1 res;
    
    //recover the component a of the RLWE sample
    FFTbackward(res[0], ct[0]);
    //recover the component b of the RLWE sample
    FFTbackward(res[1], ct[1]);
    //extract e+u from b
    for (int k = 0; k < N2; k++)
      ai[k] = ((double complex) ct[0][k]) * ((double complex) sk_FFT[k]);
    FFTbackward(asi, ai);
    // Ring_ModQ sk;
    // FFTbackward(sk, sk_FFT);
    // AntiCyclicConv(asi, res[0], sk);

    for (int k = 0; k < N; k++)
      res[1][k] =  res[1][k] - asi[k];
    //eliminate the noise e from u+e where u=Q/8
    for (int k = 0; k < N; k++){ 
      res[1][k] = res[1][k] + Q/16; // mod Q
      res[1][k] = 8*(uint64_t)res[1][k]/Q;
      if (res[1][k] == 1)
          return k*q/(2*N);
      else if (res[1][k] == 7)  // Q-1
          return q/2 + k*q/(2*N);
    }
    return -1; //exceptional return 
  }

  void KeyGen(EvalKey* EK, const LWE::SecretKey LWEsk) {
    LWE::SecretKeyN FHEWsk;
    LWE::KeyGenN(FHEWsk);
    LWE::SwitchingKeyGen(EK->KSkey, LWEsk, FHEWsk);

    Ring_FFT FHEWskFFT;
    FFTforward(FHEWskFFT, FHEWsk);
    for (int i = 0; i < n; ++i)
      for (int j = 1; j < BS_base; ++j)
        for (int k = 0; k < BS_exp; ++k) {
          EK->BSkey[i][j][k] = (ct_FFT*) fftw_malloc(sizeof(ct_FFT));
          FHEWencrypt( (*EK->BSkey[i][j][k]), FHEWskFFT, LWEsk[i] * j * BS_table[k] );
        }
  }

  void KeyGen_verbose(EvalKey* EK, Ring_FFT FHEWskFFT, const LWE::SecretKey LWEsk) { // compared with KeyGen(), it inputs FHWEsk (z) in FFT format
    LWE::SecretKeyN FHEWsk;
    FFTbackward(FHEWsk, FHEWskFFT);
    LWE::SwitchingKeyGen(EK->KSkey, LWEsk, FHEWsk); // KSkey in LWE sample

    for (int i = 0; i < n; ++i)
      for (int j = 1; j < BS_base; ++j)
        for (int k = 0; k < BS_exp; ++k) {
          EK->BSkey[i][j][k] = (ct_FFT*) fftw_malloc(sizeof(ct_FFT)); // BSkey in FFT domain
          FHEWencrypt( (*EK->BSkey[i][j][k]), FHEWskFFT, LWEsk[i] * j * BS_table[k] );
          // int m = FHEW::FHEWdecrypt((*EK->BSkey[i][j][k])[1], FHEWskFFT);
          // if ( ((LWEsk[i]*j*BS_table[k])%q+q)%q != m ){
          //   cout << "BSkey generate errors\n";
          //   cout << "i = " << i << " j = " << j << " k = " << k << "\n";
          // }
        }
  }

  void fwrite_ek(const EvalKey& EK, FILE* f) {
    // Write bootstrapping key
    for (int i = 0; i < n; ++i)      
      for (int j = 1; j < BS_base; ++j)
        for (int k = 0; k < BS_exp; ++k) 
          assert(fwrite(EK.BSkey[i][j][k], sizeof(ct_FFT), 1, f));
    // Write switching key
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < KS_base; ++j)
        for (int k = 0; k < KS_exp; ++k)
          assert(fwrite(EK.KSkey[i][j][k], sizeof(LWE::CipherTextQ), 1, f));
  }

  EvalKey* fread_ek(FILE* f) {
    EvalKey* EK = new EvalKey;
    // Read bootstrapping key
    for (int i = 0; i < n; ++i)
      for (int j = 1; j < BS_base; ++j)
        for (int k = 0; k < BS_exp; ++k) {
          EK->BSkey[i][j][k] = (ct_FFT*) fftw_malloc(sizeof(ct_FFT));
          assert(fread(EK->BSkey[i][j][k], sizeof(ct_FFT), 1, f));
        }
    // Read switching key
    cerr << "BSKey Read. \n";
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < KS_base; ++j)
        for (int k = 0; k < KS_exp; ++k) {
          EK->KSkey[i][j][k] = new LWE::CipherTextQ;
          assert(fread(EK->KSkey[i][j][k], sizeof(LWE::CipherTextQ), 1, f));
    }
    return EK;
  }

  
  
  void AddToACC(ct_FFT1 ACC, ct_FFT C) {
    ct_ModQ1 ct;
    dct_ModQ1 dct;
    dct_FFT1 dctFFT;
    // Decompose_ct(dctFFT, ACC);
    for (int j = 0; j < 2; ++j)
      FFTbackward(ct[j], ACC[j]);
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < N; ++k) {
        ZmodQ t = ct[j][k] * v_inverse;
        // int sign = t < 0 ? -1 : 1;
        // int t_abs = t * sign;
        for (int l = 0; l < K; ++l) {
          ZmodQ r = (t << g_bits_32[l]) >> g_bits_32[l];
               // if ((k ==0) && (j==0)){
               //  cout << "r = " << r << ",";}
          t = (t-r) >> g_bits[l];
          dct[j+2*l][k] = r;


          // ZmodQ r = (t_abs << g_bits_32[l]) >> g_bits_32[l];
          // if ((l==2) && (k ==0 ))
          //   cout << r << ",";
          // t_abs = (t_abs-r) >> g_bits[l];
          // dct[j+2*l][k] = r*sign;
        }
      }
    for (int j = 0; j < K2; ++j)
      FFTforward(dctFFT[j], dct[j]);
    // Mult_dct_ct(ACC, dct, C);
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < N2; ++k) {
        ACC[j][k] = (double complex) 0.0;
        for (int l = 0; l < K2; ++l)
          ACC[j][k] += ((double complex) dctFFT[l][k]) * ((double complex) C[l][j][k]);
      }
    // Use Convolution directly
    // Ring_ModQ tmp;
    // ct_ModQ1 ACC1;
    // ct_ModQ C1;
    // for (int j = 0; j < 2; j++)
    //   for (int k = 0; k < N; k++)
    //     ACC1[j][k] = 0;
    // for (int j = 0; j < 2; j++)
    //   for (int l = 0; l < K2; l++){
    //     FFTbackward(C1[l][j], C[l][j]);
    //     AntiCyclicConv(tmp, dct[l], C1[l][j]);
    //     for (int k = 0; k < N; k++)
    //       ACC1[j][k] += tmp[k];
    //   }
    // for (int j = 0; j < 2; j++)
    //   FFTforward(ACC[j], ACC1[j]);

  }
  
  void InitializeACC(ct_FFT1 ACC, int m) {  // Set a ciphertext to X^m * G (encryption of m without errors)
    ct_ModQ1 res;
    int mm = (((m % q) + q) % q) * (2*N/q);             // Reduce mod q (dealing with negative number as well)
    int sign = 1;
    if (mm >= N) { mm -= N; sign = -1; }
    
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < N; ++k)
        res[j][k]=0;
    res[1][mm] += sign*vgprime[0]; // [a,as+e] + X^m *G
    for (int j = 0; j < 2; ++j)
      FFTforward(ACC[j], res[j]);
  }
  
  LWE::CipherTextQN* MemberTest(Ring_FFT t, ct_FFT1 C) {
    Ring_FFT temp;
    Ring_ModQ temp_ModQ;
    
    LWE::CipherTextQN* ct = new LWE::CipherTextQN;
    // for (int i = 0; i < N2; ++i)
    //   temp[i] = conj(((double complex) C[0][i]) * ((double complex)t[i]));  // Compute t*a
    // FFTbackward(temp_ModQ, temp);
    Ring_ModQ t1, C0, C1;
    for (int i = 0; i < N2; ++i)
      temp[i] = conj((double complex) C[0][i]);  // Compute a(1/x)
    FFTbackward(C0, temp);
    for (int i = 0; i < N2; ++i)
      temp[i] = conj((double complex) t[i]);  // Compute t(1/x)
    FFTbackward(t1, temp);    
    AntiCyclicConv(temp_ModQ, C0, t1); // Compute a(1/x)*t(1/x)

    for (int i = 0; i < N; ++i) ct->a[i] = temp_ModQ[i];

    // for (int i = 0; i < N2; ++i)
    //   temp[i] = ((double complex) C[1][i]) * ((double complex) t[i]);
    // FFTbackward(temp_ModQ, temp);
    FFTbackward(C1, C[1]);
    FFTbackward(t1, t);
    AntiCyclicConv(temp_ModQ, C1, t1);
    ct->b = v+temp_ModQ[0]; 
    return ct;    
  }


  
  void HomGate(LWE::CipherText* res, const BinGate gate, const EvalKey& EK, const LWE::CipherText& ct1, const LWE::CipherText& ct2) {
    LWE::CipherText e12;

    for (int i = 0; i < n; ++i){
      if (((ct1.a[i] - ct2.a[i]) %q) && ((ct1.a[i] + ct2.a[i]) %q))
        // Ciphertexts are neither equal neither a NOT of each other
        break;
      if (i==n-1){
        cerr << "ERROR: Please only use independant ciphertexts as inputs." << endl;
        exit(1);
      }
    }

    for (int i = 0; i < n; ++i)
      e12.a[i] = (2*q - (ct1.a[i] + ct2.a[i])) % q;

    e12.b  =  GateConst[gate] - (ct1.b + ct2.b) % q;

    ct_FFT1 ACC;
    InitializeACC(ACC, (e12.b + q/4) % q);
    for (int i = 0; i < n; ++i) {
      int a = (q - e12.a[i] % q) % q;
      for (int k = 0; k < BS_exp; ++k, a /= BS_base) {
        int a0 = a % BS_base;
        if (a0) AddToACC(ACC, *(EK.BSkey[i][a0][k]));
      }
    }
    LWE::CipherTextQN *eQN = MemberTest(t_TestMSB, ACC);
    LWE::CipherTextQ eQ;
    LWE::KeySwitch(&eQ, EK.KSkey, *eQN);
    LWE::ModSwitch(res, eQ);
    delete eQN;
  }

  void HomNAND(LWE::CipherText* res, const EvalKey& EK, const LWE::CipherText& ct1, const LWE::CipherText& ct2) {
    HomGate(res, NAND, EK, ct1, ct2);
  }

  void HomNOT(LWE::CipherText* res, const LWE::CipherText& ct) {
    for (int i = 0; i < n; ++i)
      res->a[i] = (q - ct.a[i]) % q;
    res->b  =  (-ct.b + 5*q/4) % q;
  }

}
