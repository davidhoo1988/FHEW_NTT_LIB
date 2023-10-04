#ifndef FHEW_H
#define FHEW_H

#include "params.h"
#include "LWE.h"

namespace FHEW {
  extern Ring_NTT t_TestMSB;

  typedef Ring_NTT  ct_NTT[K2][2];  // Ciphertext in NTT form
  typedef ct_NTT* BootstrappingKey[n][BS_base][BS_exp];
  typedef struct {
    BootstrappingKey BSkey;
    LWE::SwitchingKey KSkey;
  } EvalKey;
  typedef Ring_NTT ct_NTT1[2]; // C[1] entry of ct_NTT ciphertext
  typedef Ring_ModQ ct_ModQ[K2][2];   // Ciphertext in coefficient form
  typedef Ring_ModQ ct_ModQ1[2];   // Entry C[1] of ct_ModQ
  typedef Ring_ModQ dct_ModQ[K2][K2]; // Decomposed Ciphertext in coeff form
  typedef Ring_ModQ dct_ModQ1[K2]; // Entry C[1] of dct_ModQ
  typedef Ring_NTT  dct_NTT[K2][K2];  // Decomposed Ciphertext in NTT form
  typedef Ring_NTT  dct_NTT1[K2];  // Entry C[1] of dct_NTT

  void Setup();
  void FHEWencrypt(ct_NTT ct, Ring_NTT sk_NTT, int m);
  int FHEWdecrypt(ct_NTT1 ct, Ring_NTT sk_NTT);
  void KeyGen(EvalKey* EK, const LWE::SecretKey LWEsk);
  void KeyGen_verbose(EvalKey* EK, Ring_NTT FHEWskFFT, const LWE::SecretKey LWEsk); // For debugging only, compared with KeyGen(), it inputs FHWEsk (z) in NTT format
  void HomGate(LWE::CipherText* res, const BinGate gate, const EvalKey& EK, const LWE::CipherText& ct1, const LWE::CipherText& ct2);
  void HomNAND(LWE::CipherText* res, const EvalKey& EK, const LWE::CipherText& ct1, const LWE::CipherText& ct2);
  void HomNOT(LWE::CipherText* res, const LWE::CipherText& ct);
  void fwrite_ek(const EvalKey& EK, FILE* f);
  EvalKey* fread_ek(FILE* f);

  //some internal functions, can be hidden
  void AddToACC(ct_NTT1 ACC, ct_NTT C);
  void InitializeACC(ct_NTT1 ACC, int m);
  LWE::CipherTextQN* MemberTest(Ring_NTT t, ct_NTT1 C);
}

#endif
