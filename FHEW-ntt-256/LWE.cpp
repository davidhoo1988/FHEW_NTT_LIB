#include <cstdlib>
#include <iostream>
#include "LWE.h"
#include "params.h"
#include <cassert>

using namespace std;
namespace LWE {

  void KeyGen(SecretKey sk) {
    KeyGenRestart:
    int s=0, ss=0;
    for (int i = 0; i < n; ++i) {

          sk[i] = Sample(Chi_Binary);
          s+= sk[i];
          ss+= abs(sk[i]);
        }
      if (abs(s)>5) goto KeyGenRestart;
      if (abs(ss - n/2)>5) goto KeyGenRestart;
  }

  void KeyGenN(SecretKeyN sk) {
    int tp = 0;
    for (int i = 0; i < N; ++i){ 
      tp = Sample(Chi1);  
      sk[i] = tp;
      //sk[i] = tp < 0 ? Q+tp : tp; //make every entry in sk small {0,1}, -1 is problematic for key-switching procedure!
    }
  }
  
  void Encrypt(CipherText* ct, const SecretKey sk, int m) {
    ct->b = (m % 4) * q / 4;
    ct->b += Sample(Chi3);
    for (int i = 0; i < n; ++i)	{
      ct->a[i] = rand() % q;
      ct->b = (ct->b + ct->a[i] * sk[i]) % q;
    }
  }

  void Encrypt2(CipherText* ct, const SecretKey sk, int m) { // LWE for 2/q and s
    ct->b = (m % 2) * q / 2;
    ct->b += Sample(Chi3);
    for (int i = 0; i < n; ++i) {
      ct->a[i] = rand() % q;
      ct->b = (ct->b + ct->a[i] * sk[i]) % q;
    }
  }

  void EncryptQ(CipherTextQ* ct, const SecretKey sk, int m) {
    ct->b = (m % 4) * (Q / 4 + Q % 4);
    // cerr << "normalized message: " << ct->b << endl;
    ct->b = ModAddQ(ct->b, Sample(Chi3));
    for (int i = 0; i < n; ++i) {
      ct->a[i] = ModMulQ(rand(), 1) ;
      ct->b = ModAddQ(ct->b, ModMulQ(ct->a[i], sk[i]));
    }
  }

  void EncryptQN(CipherTextQN* ct, const SecretKeyN sk, int m) {
    ct->b = (m % 4) * (Q / 4 + Q % 4);
    // cerr << "normalized message: " << ct->b << endl;
    ct->b = ModAddQ(ct->b, Sample(Chi3));
    for (int i = 0; i < N; ++i) {
      ct->a[i] = ModMulQ(rand(), 1) ;
      ct->b = ModAddQ(ct->b, ModMulQ(ct->a[i], sk[i]));
    }
  }

  int Decrypt(const SecretKey sk, const CipherText& ct) {
    int r = ct.b;
    for (int i = 0; i < n; ++i) r -= ct.a[i] * sk[i];
    r = ((r % q) + q + q/8) % q;
    return 4*r/q;    
  }

  int Decrypt2(const SecretKey sk, const CipherText& ct) { // LWE for 2/q and s
    int r = ct.b;
    for (int i = 0; i < n; ++i) r -= ct.a[i] * sk[i];
    r = ((r % q) + q + q/4) % q;
    return 2*r/q;    
  }

  int DecryptQ(const SecretKey sk, const CipherTextQ& ct) {
    ZmodQ r = ct.b;
    for (int i = 0; i < n; ++i) r = ModSubQ(r, ModMulQ(ct.a[i], sk[i]));
    // cerr << "normalized message1: " << r << endl;
    r = floor(.5 + ((double)4.0/Q)*r);
    return (r+4)%4;    
  }

  int DecryptQN(const SecretKeyN sk, const CipherTextQN& ct) {
    ZmodQ r = ct.b;
    for (int i = 0; i < N; ++i) r = ModSubQ(r, ModMulQ(ct.a[i], sk[i]));
    // cerr << "normalized message1: " << r << endl;
    r = floor(.5 + ((double)4.0/Q)*r); 
    return (r+4)%4;     
  }

  void DecryptDetail(const SecretKey sk, const CipherText& ct) {
    int r = ct.b;
    for (int i = 0; i < n; ++i) r -= ct.a[i] * sk[i];
    r = ((r % q) + q + q/8) % q;
    int m = 4*r/q;
    cerr << "\t Value " << r - q/8 << "\t Decoded as " << m << " * " << q/4 << "\t + " << r - m * q/4 - q/8<< endl;
    cout << r - m * q/4 - q/8 << ", "; 
  }

    
  int round_qQ(ZmodQ v) {
    return floor(.5 + (double) v * (double) q / (double) Q) + q % q;
  }
  
  void ModSwitch(CipherText* ct, const CipherTextQ& c) {
    for (int i = 0; i < n; ++i) 
      ct->a[i] = round_qQ(c.a[i]);	
    ct->b = round_qQ(c.b);
  }

  void SwitchingKeyGen(SwitchingKey res, const SecretKey new_sk, const SecretKeyN old_sk) {
    //int old_sk_entry = 0;
    for (int i = 0; i < N; ++i) 
      for (int j = 0; j < KS_base; ++j)
  	for (int k = 0; k < KS_exp; ++k) {
  	  CipherTextQ* ct = new CipherTextQ;
      // if (old_sk[i] > Q/2) //transform old_sk back to [-q/2, q/2]
      //    old_sk_entry = old_sk[i]-Q;   
      // else 
      //    old_sk_entry = old_sk[i];
  	  ct->b = ModAddQ(ModMulQ(KS_table[k], ModMulQ(-old_sk[i], j)), Sample(Chi1)); 
  	  for (int l = 0; l < n; ++l) {
  	    ct->a[l] = ModMulQ(rand(),1);
        ct->b = ModAddQ(ct->b, ModMulQ(ct->a[l], new_sk[l]));
  	  }
  	  res[i][j][k] = ct;
  	}
  }
  
  void KeySwitch(CipherTextQ* res, const SwitchingKey K, const CipherTextQN& ct) {
    for (int k = 0; k < n; ++k) res->a[k] = 0;
    res->b = ct.b;
    for (int i = 0; i < N; ++i) {
      uZmodQ a = 0;
      if (ct.a[i] <= 0)
        a = -ct.a[i];
      else
        a = Q-ct.a[i];  
      for (int j = 0; j < KS_exp; ++j, a /= KS_base) {
      	uZmodQ a0 = a % KS_base;
      	for (int k = 0; k < n; ++k)
          res->a[k] = ModSubQ(res->a[k], (K[i][a0][j])->a[k]);
        res->b = ModSubQ(res->b, (K[i][a0][j])->b);
      }
    } 
  }

  ZmodQ ModMulQ(ZmodQ x, ZmodQ y){
    int64_t xy = int64_t(x) * int64_t(y);
    // while (xy > Q){
    //   uint64_t xy0 = xy & 0xffffffff;
    //   uint64_t xy1 = (xy>>32) & 0xffffffff;  
    //   xy = xy0 + 10239*xy1; // recall Q = 2^32 - 10239 
    //   //cout << "xy = " << xy << "\n";
    // }
    xy = xy % Q; 
    if (xy > int64_t(Q)/2)
      xy = xy - Q;
    else if (xy < -(int64_t(Q)/2)) 
      xy = xy + Q;  
    return ZmodQ(xy);
  }

  // ZmodQ ModMulQ(ZmodQ x, uZmodQ y){
  //   int64_t xy = int64_t(x) * int64_t(y);
  //   // uint64_t xy0 = xy & 0xffffffff;
  //   // uint64_t xy1 = (xy>>32) & 0xffffffff;  
  //   // xy = xy0 + 10239*xy1; // recall Q = 2^32 - 10239
  //   xy = xy % Q;
  //   return ZmodQ(xy);
  // }

  ZmodQ ModAddQ(ZmodQ x, ZmodQ y){
    int64_t xy = int64_t(x) + int64_t(y);
    //xy = xy % Q;
    if (xy > int64_t(Q)/2)
      xy = xy - Q;
    else if (xy < -(int64_t(Q)/2)) 
      xy = xy + Q; 
    return ZmodQ(xy);
  }

  // ZmodQ ModAddQ(ZmodQ x, uZmodQ y){
  //   int64_t xy = int64_t(x) + int64_t(y);
  //   if (xy > Q/2)
  //     xy = xy - Q;
  //   else if (xy < -int64_t(Q)/2+1) 
  //     xy = xy + Q; 
  //   return ZmodQ(xy);
  // }

  ZmodQ ModSubQ(ZmodQ x, ZmodQ y){
    int64_t xy = int64_t(x) - int64_t(y);
    //xy = xy % Q;
    if (xy > int64_t(Q)/2)
      xy = xy - Q;
    else if (xy < -(int64_t(Q)/2)) 
      xy = xy + Q; 
    return ZmodQ(xy);
  }

}
