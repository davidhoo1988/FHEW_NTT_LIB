#ifndef FHEW_H
#define FHEW_H

#include "params.h"
#include "LWE.h"

namespace FHEW {
  extern Ring_FFT t_TestMSB;

  typedef Ring_FFT  ct_FFT[K2][2];  // Ciphertext in FFT form
  typedef ct_FFT* BootstrappingKey[n][BS_base][BS_exp];
  typedef struct {
    BootstrappingKey BSkey;
    LWE::SwitchingKey KSkey;
  } EvalKey;
  typedef Ring_FFT ct_FFT1[2]; // C[1] entry of ct_FFT ciphertext
  typedef Ring_ModQ ct_ModQ[K2][2];   // Ciphertext in coefficient form
  typedef Ring_ModQ ct_ModQ1[2];   // Entry C[1] of ct_ModQ
  typedef Ring_ModQ dct_ModQ[K2][K2]; // Decomposed Ciphertext in coeff form
  typedef Ring_ModQ dct_ModQ1[K2]; // Entry C[1] of dct_ModQ
  typedef Ring_FFT  dct_FFT[K2][K2];  // Decomposed Ciphertext in FFT form
  typedef Ring_FFT  dct_FFT1[K2];  // Entry C[1] of dct_FFT

  void Setup();
  void FHEWencrypt(ct_FFT ct, Ring_FFT sk_FFT, int m);
  int FHEWdecrypt(ct_FFT1 ct, Ring_FFT sk_FFT);
  void KeyGen(EvalKey* EK, const LWE::SecretKey LWEsk);
  void KeyGen_verbose(EvalKey* EK, Ring_FFT FHEWskFFT, const LWE::SecretKey LWEsk); // For debugging only, compared with KeyGen(), it inputs FHWEsk (z) in FFT format
  void HomGate(LWE::CipherText* res, const BinGate gate, const EvalKey& EK, const LWE::CipherText& ct1, const LWE::CipherText& ct2);
  void HomNAND(LWE::CipherText* res, const EvalKey& EK, const LWE::CipherText& ct1, const LWE::CipherText& ct2);
  void HomNOT(LWE::CipherText* res, const LWE::CipherText& ct);
  void fwrite_ek(const EvalKey& EK, FILE* f);
  EvalKey* fread_ek(FILE* f);

  //some internal functions, can be hidden
  void AddToACC(ct_FFT1 ACC, ct_FFT C);
  void InitializeACC(ct_FFT1 ACC, int m);
  LWE::CipherTextQN* MemberTest(Ring_FFT t, ct_FFT1 C);
}

#endif
