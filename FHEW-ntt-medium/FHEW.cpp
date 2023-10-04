#include <iostream>
#include <cstdlib>
#include "FHEW.h"
#include "FFT.h"
#include <cassert>

using namespace std;
namespace FHEW {  
  Ring_NTT t_TestMSB;

  void Setup() {
    NTTsetup();
    Ring_ModQ tmsb;
    tmsb[0]=-1; // -1 in F_Q
    for (int i = 1; i < N; ++i) tmsb[i]=1;
    NTTforward(t_TestMSB, tmsb);
  }
  
  void FHEWencrypt(ct_NTT ct, Ring_NTT sk_NTT, int m) { // Used by BootstrappingKeyGen
    Ring_NTT ai;
    ct_ModQ res;
    int mm = (((m % q) + q) % q) * (2*N/q);             // Reduce mod q (dealing with negative number as well)
    int sign = 1;
    if (mm >= N) { mm -= N; sign = -1; }
    
    for (int i = 0; i < K2; ++i) {
      for (int k = 0; k < N; ++k) res[i][0][k] = LWE::ModMulQ(rand(), 1); // % Q  LWE::ModMulQ(rand(), 1)
      NTTforward(ai, res[i][0]);
      for (int k = 0; k < N; ++k) 
        ai[k] = LWE::ModMulQ(ai[k], sk_NTT[k]);
      NTTbackward(res[i][1], ai);
      for (int k = 0; k < N; ++k) res[i][1][k] = LWE::ModAddQ(res[i][1][k], Sample(Chi1));    // Add error [a,as+e]
    }
    // // Use Convolution directly
    // Ring_ModQ sk;
    // NTTbackward(sk, sk_NTT);
    // for (int i = 0; i < K2; ++i) {
    //   for (int k = 0; k < N; ++k) res[i][0][k] = (-1*rand())%(1<<16);  // % Q
    //   AntiCyclicConv(res[i][1], res[i][0], sk);
    //   for (int k = 0; k < N; ++k) res[i][1][k] += 0;    // Add error [a,as+e]
    // }
    
    for (int i = 0; i < K; ++i) {
      res[2*i  ][0][mm] = LWE::ModAddQ(res[2*i  ][0][mm], LWE::ModMulQ(vgprime[i], sign)); // Add G Multiple
      res[2*i+1][1][mm] = LWE::ModAddQ(res[2*i+1][1][mm], LWE::ModMulQ(vgprime[i], sign)); // [a,as+e] + X^m *G
    }
    for (int i = 0; i < K2; ++i)
      for (int j = 0; j < 2; ++j)
        NTTforward(ct[i][j], res[i][j]);
  }
  
  int FHEWdecrypt(ct_NTT1 ct, Ring_NTT sk_NTT){
    Ring_NTT ai;
    Ring_ModQ asi;
    ct_ModQ1 res;
    
    //recover the component a of the RLWE sample
    NTTbackward(res[0], ct[0]);
    //recover the component b of the RLWE sample
    NTTbackward(res[1], ct[1]);
    //extract e+u from b
    for (int k = 0; k < N; k++)
      ai[k] = LWE::ModMulQ(ct[0][k], sk_NTT[k]);
    NTTbackward(asi, ai);
    // Ring_ModQ sk;
    // NTTbackward(sk, sk_NTT);
    // AntiCyclicConv(asi, res[0], sk);
    for (int k = 0; k < N; k++)
      res[1][k] = LWE::ModSubQ(res[1][k], asi[k]);
    // cout << "ACC after decrytion: " << "\n";
    // for (int k = 0; k < N; k++)
    //   cout << res[1][k] << "\t";
    //eliminate the noise e from u+e where u \approx Q/8
    for (int k = 0; k < N; k++){ 
      //res[1][k] = LWE::ModAddQ(res[1][k], v/2); // + Q/16 mod Q
      //res[1][k] = 8*(uint64_t)res[1][k]/Q;
      res[1][k] = floor(.5 + double(8)*res[1][k]/double(Q));
      if (res[1][k] == 1)
          return k*q/(2*N);
      else if (res[1][k] == -1)  // -1
          return q/2 + k*q/(2*N);
    }
    return -1; //exceptional return 
  }

  void KeyGen(EvalKey* EK, const LWE::SecretKey LWEsk) {
    LWE::SecretKeyN FHEWsk;
    LWE::KeyGenN(FHEWsk);
    LWE::SwitchingKeyGen(EK->KSkey, LWEsk, FHEWsk);

    Ring_NTT FHEWskNTT;
    NTTforward(FHEWskNTT, FHEWsk);
    for (int i = 0; i < n; ++i)
      for (int j = 1; j < BS_base; ++j)
        for (int k = 0; k < BS_exp; ++k) {
          EK->BSkey[i][j][k] = (ct_NTT*) malloc(sizeof(ct_NTT));
          FHEWencrypt( (*EK->BSkey[i][j][k]), FHEWskNTT, LWEsk[i] * j * BS_table[k] );
        }
  }

  void KeyGen_verbose(EvalKey* EK, Ring_NTT FHEWskNTT, const LWE::SecretKey LWEsk) { // compared with KeyGen(), it inputs FHWEsk (z) in NTT format
    LWE::SecretKeyN FHEWsk;
    NTTbackward(FHEWsk, FHEWskNTT);
    LWE::SwitchingKeyGen(EK->KSkey, LWEsk, FHEWsk); // KSkey in LWE sample

    for (int i = 0; i < n; ++i)
      for (int j = 1; j < BS_base; ++j)
        for (int k = 0; k < BS_exp; ++k) {
          EK->BSkey[i][j][k] = (ct_NTT*) malloc(sizeof(ct_NTT)); // BSkey in NTT domain
          FHEWencrypt( (*EK->BSkey[i][j][k]), FHEWskNTT, LWEsk[i] * j * BS_table[k] );
          // int m = FHEW::FHEWdecrypt((*EK->BSkey[i][j][k])[1], FHEWskNTT);
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
          assert(fwrite(EK.BSkey[i][j][k], sizeof(ct_NTT), 1, f));
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
          EK->BSkey[i][j][k] = (ct_NTT*) malloc(sizeof(ct_NTT));
          assert(fread(EK->BSkey[i][j][k], sizeof(ct_NTT), 1, f));
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

  
  
  void AddToACC(ct_NTT1 ACC, ct_NTT C) {
    ct_ModQ1 ct;
    dct_ModQ1 dct;
    dct_NTT1 dctNTT;
    // Decompose_ct(dctNTT, ACC);
    for (int j = 0; j < 2; ++j)
      NTTbackward(ct[j], ACC[j]);
    for (int j = 0; j < 2; ++j){
      for (int k = 0; k < N; ++k) {
        ZmodQ t = LWE::ModMulQ(ct[j][k], v_inverse);
        //int sign = t < 0 ? -1 : 1;
        //int t_abs = t * sign;
        for (int l = 0; l < K; ++l) {
          ZmodQ r = (t << g_bits_31[l]) >> g_bits_31[l];
               // if ((k ==0) && (j==0)){
               //  cout << "r = " << r << ",";}
          t = (t-r) >> g_bits[l];
          dct[j+2*l][k] = r;
          //ZmodQ r = (t_abs << g_bits_32[l]) >> g_bits_32[l];
          // if ((l==2) && (k ==0 ))
          //   cout << r << ",";
          //t_abs = (t_abs-r) >> g_bits[l];
          //dct[j+2*l][k] = r*sign;
        }
      }
    }
      //check the correctness of decompose()
      for (int j = 0; j < 2; j++)
        for (int k = 0; k < N; k++){
          ZmodQ tmp = 0;
          for (int l = 0; l < K; l++)
            tmp = LWE::ModAddQ(tmp, LWE::ModMulQ( dct[j+2*l][k], 1<<(l*g_bits[0]) ));
            //tmp += dct[j+2*l][k] * vgprime[l];
          if (tmp != LWE::ModMulQ(ct[j][k], v_inverse))
            cout << "error: " << "j,k = " << j << ", " << k << "\t" << LWE::ModMulQ(ct[j][k], v_inverse) << " <---> " << tmp << "\n";
        }

    for (int j = 0; j < K2; ++j)
      NTTforward(dctNTT[j], dct[j]);
    // Mult_dct_ct(ACC, dct, C);
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < N; ++k) {
        ACC[j][k] = (ZmodQ) 0;
        for (int l = 0; l < K2; ++l)
          ACC[j][k] = LWE::ModAddQ(ACC[j][k], LWE::ModMulQ(dctNTT[l][k], C[l][j][k]));
      }
    // // Use Convolution directly
    // Ring_ModQ tmp;
    // ct_ModQ1 ACC1;
    // ct_ModQ C1;
    // for (int j = 0; j < 2; j++)
    //   for (int k = 0; k < N; k++)
    //     ACC1[j][k] = 0;
    // for (int j = 0; j < 2; j++)
    //   for (int l = 0; l < K2; l++){
    //     NTTbackward(C1[l][j], C[l][j]);
    //     AntiCyclicConv(tmp, dct[l], C1[l][j]);
    //     for (int k = 0; k < N; k++)
    //       ACC1[j][k] = LWE::ModAddQ(ACC1[j][k], tmp[k]);
    //   }
    // for (int j = 0; j < 2; j++)
    //   NTTforward(ACC[j], ACC1[j]);

  }
  
  void InitializeACC(ct_NTT1 ACC, int m) {  // Set a ciphertext to X^m * G (encryption of m without errors)
    ct_ModQ1 res;
    int mm = (((m % q) + q) % q) * (2*N/q);             // Reduce mod q (dealing with negative number as well)
    int sign = 1;
    if (mm >= N) { mm -= N; sign = -1; }
    
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < N; ++k)
        res[j][k] = 0;
    res[1][mm] = LWE::ModAddQ(res[1][mm], LWE::ModMulQ(vgprime[0], sign)); // [a,as+e] + X^m *G
    for (int j = 0; j < 2; ++j)
      NTTforward(ACC[j], res[j]);
  }

  void InitializeACC2(ct_NTT1 ACC, int m) {  // Set a ciphertext to testMSB * X^m * G (encryption of m without errors)
    ct_ModQ1 res;

    int mm = (((m % q) + q) % q) * (2*N/q);             // Reduce mod q (dealing with negative number as well)
    int sign = 1;
    if (mm >= N) { mm -= N; sign = -1; }
    
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < N; ++k)
        res[j][k] = 0;
    res[1][mm] = LWE::ModAddQ(res[1][mm], LWE::ModMulQ(vgprime[0], sign)); // [0,0] + X^m *G
    for (int j = 0; j < 2; ++j)
      NTTforward(ACC[j], res[j]);

    for (int i = 0; i < N; ++i){ // Compute t*a
      ACC[0][i] = LWE::ModMulQ(ACC[0][i], t_TestMSB[i]);
      ACC[1][i] = LWE::ModMulQ(ACC[1][i], t_TestMSB[i]);  
    }
  }

  LWE::CipherTextQN* MemberTest(Ring_NTT t, ct_NTT1 C) {
    Ring_NTT temp;
    Ring_ModQ temp_ModQ;
    
    LWE::CipherTextQN* ct = new LWE::CipherTextQN;
    for (int i = 0; i < N; ++i)
      temp[i] = LWE::ModMulQ(C[0][i], t[i]);  // Compute t*a
    NTTbackward(temp_ModQ, temp);   
    for (int i = 0; i < N; ++i)  // Compute Const(t*a)
      if (i == 0)
        ct->a[i] = temp_ModQ[i];
      else
        ct->a[i] = LWE::ModMulQ(temp_ModQ[N-i], -1);

    // Ring_ModQ t1, C1;
    // NTTbackward(C1, C[1]);
    // NTTbackward(t1, t);
    // AntiCyclicConv(temp_ModQ, C1, t1); // Compute t*b
    for (int i = 0; i < N; ++i)
      temp[i] = LWE::ModMulQ(C[1][i], t[i]);  // Compute t*a
    NTTbackward(temp_ModQ, temp);  

    ct->b = LWE::ModAddQ(v, temp_ModQ[0]); 
    return ct;    
  }

  LWE::CipherTextQN* MemberTest2(ct_NTT1 C) { //extract constant term from the RLWE encryption as an LWE sample
    Ring_ModQ temp_ModQ;
    
    LWE::CipherTextQN* ct = new LWE::CipherTextQN;
    NTTbackward(temp_ModQ, C[0]);   
    for (int i = 0; i < N; ++i)  // Compute Const(a)
      if (i == 0)
        ct->a[i] = temp_ModQ[i];
      else
        ct->a[i] = LWE::ModMulQ(temp_ModQ[N-i], -1);

    // Ring_ModQ t1, C1;
    // NTTbackward(C1, C[1]);
    // NTTbackward(t1, t);
    // AntiCyclicConv(temp_ModQ, C1, t1); // Compute t*b
    NTTbackward(temp_ModQ, C[1]);  

    ct->b = LWE::ModAddQ(v, temp_ModQ[0]); 
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

    ct_NTT1 ACC;
    //InitializeACC(ACC, (e12.b + q/4) % q);
    InitializeACC2(ACC, (e12.b + q/4) % q);
    for (int i = 0; i < n; ++i) {
      int a = (q - e12.a[i] % q) % q;
      for (int k = 0; k < BS_exp; ++k, a /= BS_base) {
        int a0 = a % BS_base;
        if (a0) AddToACC(ACC, *(EK.BSkey[i][a0][k]));
      }
    }
    //LWE::CipherTextQN *eQN = MemberTest(t_TestMSB, ACC);
    LWE::CipherTextQN *eQN = MemberTest2(ACC);
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
