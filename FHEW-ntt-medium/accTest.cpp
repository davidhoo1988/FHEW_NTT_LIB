#include <iostream>
#include <cstdlib>
#include "LWE.h"
#include "FHEW.h"
#include "FFT.h"
#include "distrib.h"


using namespace std;

LWE::SecretKeyN FHEWsk_z;
Ring_NTT FHEWsk_zNTT;
FHEW::ct_NTT FHEWctNTT;
FHEW::ct_NTT1* FHEWctNTT1;

void help(char* cmd) {
  cerr << "\nusage: " << cmd << " n\n\n" 
  << "  Generate one secret keys sk_z for E_z(), and repeat the following test n times:\n"
  << "  prepare an LWE sample e12=(a,b) under sk_s and q/2"
  << "  perform ACC to get E_z(q/4+b-as)=E_z(q/4+q/2*m+e) with |e|<q/4 (decrypt E_z(.) to find 0<q/4+q/2*m+e<q/2 if m=0 or q/2<q/4+q/2*m+e<q if m=1)"
  << "  finally perform the homomorphic MSB extractor to get LWE(m) (encrypted under z and Q/4) from E_z(q/4+q/2*m+e)"
  << "\n If any of the tests fails, print ERROR and stop immediately.\n\n";
  exit(0);
}



int main(int argc, char *argv[]) {
  if (argc != 2) help(argv[0]);
  int count = atoi(argv[1]); 
  srand(time(NULL));
  int flag = 0;

  cerr << "Setting up GSW secret key \n";
  LWE::KeyGenN(FHEWsk_z);
  cerr << "Setting up FFTW initialization\n";
  FHEW::Setup();
  NTTforward(FHEWsk_zNTT, FHEWsk_z);

  cerr << "Setting up LWE secret key\n";
  LWE::SecretKey LWEsk_s;
  LWE::KeyGen(LWEsk_s);
  cerr << "Setting up ACC BSkey\n";
  FHEW::EvalKey EK;
  FHEW::KeyGen_verbose(&EK, FHEWsk_zNTT, LWEsk_s);

  cerr << "Testing ACC\n";
  for (int i = 0; i < count; i++){
    // TASk1
    int message = rand()%2; //testbench 400 (rand()%2)
    cout << message << "\n";

    FHEW::FHEWencrypt(FHEWctNTT, FHEWsk_zNTT, message);
    // int m1 = FHEW::FHEWdecrypt(FHEWctNTT[1], FHEWsk_zNTT);
    // cout << m1 << "\n";

    // TASK2
    LWE::CipherText e12; 
    FHEW::ct_NTT1 ACC; // RLWE in NTT domain
    LWE::Encrypt2(&e12, LWEsk_s, message); // e12 = LWE(message) under s and q/2

    int message1 = (e12.b+q/4)%q; //testbench 312  (e12.b+q/4)
    cout << message1 << "\n";
    FHEW::InitializeACC(ACC, message1); // ACC := FHEW(b+q/4) under z and Q (note that b+q/4 \in Z_q)
    // for (int j = 0; j < 100; j++)
    //   FHEW::AddToACC(ACC, FHEWctNTT);
    for (int i = 0; i < n; ++i) {
      int a = (q - e12.a[i] % q) % q; // decompose(-a)
      for (int k = 0; k < BS_exp; ++k, a /= BS_base) {
        int a0 = a % BS_base;
        if (a0) FHEW::AddToACC(ACC, *(EK.BSkey[i][a0][k]));
      }
    }
    //int m1 = FHEW::FHEWdecrypt(ACC, FHEWsk_zNTT);
    LWE::CipherTextQN* ct_msb = FHEW::MemberTest(FHEW::t_TestMSB, ACC);
    //int m1 = LWE::DecryptQN(FHEWsk_z, *ct_msb);
    LWE::CipherTextQ eQ;
    LWE::CipherText res;
    LWE::KeySwitch(&eQ, EK.KSkey, *ct_msb);
    LWE::ModSwitch(&res, eQ);
    int m1 = LWE::Decrypt(LWEsk_s, res);

    cout << "ACC decrypted: " << m1 << "\n"; // 0<m1<q/2 if m=0 or q/2<m1<q if m=1   
    // if ((m1 < q/2 && message == 1) || (m1 > q/2 && message == 0)){
    //   flag = -1;
    //   break;
    // }
    if ( m1 != message ){
        flag = -1;
        cout << "count = " << i << "\n";
        break;
    }
  }
  if (flag)
    cerr << "\nTests failed!\n\n";
  else
    cerr << "\nPassed all tests!\n\n";

}