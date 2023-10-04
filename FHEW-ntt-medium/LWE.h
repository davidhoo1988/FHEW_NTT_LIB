#ifndef LWE_H
#define LWE_H

#include "params.h"
#include "distrib.h"

namespace LWE {
  
  typedef struct { // n-dimensional LWE ciphertext modulo q
    int a[n]; // array of n integers modulo q 
    int b;    // integer modulo q
  } CipherText;

  typedef struct { // n-dimensional LWE ciphertext modulo Q
    ZmodQ a[n];
    ZmodQ b;
  } CipherTextQ;
  
  typedef struct { // N-dimensional LWE ciphertext modulo Q
    ZmodQ a[N];
    ZmodQ b;
  } CipherTextQN;
  
  typedef int SecretKey[n]; // n-dimensional LWE secret key 
  typedef int SecretKeyN[N]; // N-dimensional LWE secret key 
  
  void KeyGen(SecretKey sk);
  void KeyGenN(SecretKeyN sk);
  void Encrypt(CipherText* ct, const SecretKey sk, int m);
  void Encrypt2(CipherText* ct, const SecretKey sk, int m);
  void EncryptQ(CipherTextQ* ct, const SecretKey sk, int m);
  void EncryptQN(CipherTextQN* ct, const SecretKeyN sk, int m);
  int Decrypt(const SecretKey sk, const CipherText& ct);
  int Decrypt2(const SecretKey sk, const CipherText& ct);
  int DecryptQN(const SecretKeyN sk, const CipherTextQN& ct);
  int DecryptQ(const SecretKey sk, const CipherTextQ& ct);

  typedef CipherTextQ* SwitchingKey[N][KS_base][KS_exp];
  // Generate key material (SwitchingKey) required by KeySwitch to transform 
  // LWE encryptions under old_sk into LWE encryptions under new_sk
  void SwitchingKeyGen(SwitchingKey res, const SecretKey new_sk, const SecretKeyN old_sk);
  void KeySwitch(CipherTextQ* res, const SwitchingKey K, const CipherTextQN& ct);

  // Changes an LWE ciphertext modulo Q into an LWE ciphertext modulo q
  void ModSwitch(CipherText* ct, const CipherTextQ& c); 
  int round_qQ(ZmodQ v);

  // basic Mod Q arithmetic
  ZmodQ ModMulQ(ZmodQ x, ZmodQ y);

  //ZmodQ ModMulQ(ZmodQ x, uZmodQ y);
  ZmodQ ModAddQ(ZmodQ x, ZmodQ y);
  //ZmodQ ModAddQ(ZmodQ x, uZmodQ y);
  ZmodQ ModSubQ(ZmodQ x, ZmodQ y);  

// For debbugging purpose you can use the following 
//  void DecryptDetail(const SecretKey sk, const CipherText& ct); 


}




#endif
