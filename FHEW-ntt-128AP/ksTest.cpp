#include <iostream>
#include <cstdlib>
#include "LWE.h"

#include "distrib.h"

using namespace std;

LWE::SecretKeyN LWEskN;
LWE::SecretKey LWEsk;
LWE::CipherTextQN ct_old;
LWE::CipherTextQ ct_new;
LWE::CipherText ct;
LWE::SwitchingKey KSkey;

void help(char* cmd) {
  cerr << "\nusage: " << cmd << " n\n\n" 
  << "  Generate two secret keys sk_z, sk_s and keyswitching key ksk, and repeat the following test n times:\n"
  << "   - generate random bit m\n"
  << "   - compute an LWE ciphertext c1 under sk_z (mod Q)\n"
  << "   - homomorphically switching to another LWE ciphertext c2 under sk_s (mod Q)\n"
  << "   - homomorphically switching to another LWE ciphertext c3 under sk_s (mod q)\n"
  << "   - decrypt c1, c2, c3 and check correctness \n"
  << "\n If any of the tests fails, print ERROR and stop immediately.\n\n";
  exit(0);
}



int main(int argc, char *argv[]) {
  if (argc != 2) help(argv[0]);
  int count = atoi(argv[1]); 
  srand(time(NULL));
  int flag = 0;

  cerr << "Setting up LWE (old) secret key \n";
  LWE::KeyGenN(LWEskN);
  cerr << "Setting up LWE (new) secret key \n";
  LWE::KeyGen(LWEsk);
  cerr << "Setting up LWE key switching key \n";
  LWE::SwitchingKeyGen(KSkey, LWEsk, LWEskN);

  for (int i = 0; i < count; i++){
    // //Test modular arithmetic
    // int x = rand()%Q;
    // if (x > int64_t(Q)/2)
    //   x = x - Q;
    // else if (x < -(int64_t(Q)/2)+1) 
    //   x = x + Q; 

    // int y = rand()%Q;
    // if (y > int64_t(Q)/2)
    //   y = y - Q;
    // else if (y < -(int64_t(Q)/2)+1) 
    //   y = y + Q; 
    // int xy = LWE::ModMulQ(x,y);
    // if ( (xy - (int64_t(x)*int64_t(y))%Q) % Q != 0 ){
    //   cout << "i = " << i << endl;
    //   cout << "x = " << x << endl;
    //   cout << "y = " << y << endl; 
    //   cout << "xy = " << xy << endl;
    //   cout << "x*y mod Q = " << (int64_t(x)*int64_t(y))%Q << endl;
    //   cout << "fdf = " <<   xy - (int64_t(x)*int64_t(y))%Q << endl;
    //   flag = -1;
    //   break;
    // }

    // xy = LWE::ModAddQ(x,y);
    // if ( (xy - (int64_t(x)+int64_t(y))%Q) % Q != 0 ){
    //   cout << "i = " << i << endl;
    //   cout << "x = " << x << endl;
    //   cout << "y = " << y << endl; 
    //   cout << "x+y = " << xy << endl;
    //   cout << "x+y mod Q = " << (int64_t(x)+int64_t(y))%Q << endl;
    //   cout << "fdf = " <<   xy - (int64_t(x)+int64_t(y))%Q << endl;
    //   flag = -1;
    //   break;
    // }

    // xy = LWE::ModSubQ(x,y);
    // if ( (xy - (int64_t(x)-int64_t(y))%Q) % Q != 0 ){
    //   cout << "i = " << i << endl;
    //   cout << "x = " << x << endl;
    //   cout << "y = " << y << endl; 
    //   cout << "x-y = " << xy << endl;
    //   cout << "x-y mod Q = " << (int64_t(x)-int64_t(y))%Q << endl;
    //   cout << "fdf = " <<   xy - (int64_t(x)-int64_t(y))%Q << endl;
    //   flag = -1;
    //   break;
    // }

    // test LWE encryption
    int message = rand()%2;
    cout << message << "\n";

    LWE::Encrypt2(&ct, LWEsk, message);
    int m0 = LWE::Decrypt2(LWEsk, ct);
    cout << m0 << "\n";

    LWE::EncryptQN(&ct_old, LWEskN, message);
    int m = LWE::DecryptQN(LWEskN, ct_old);
    cout << m << "\n";

    LWE::KeySwitch(&ct_new, KSkey, ct_old);
    int m1 = LWE::DecryptQ(LWEsk, ct_new);
    cout << m1 << "\n";

    ModSwitch(&ct, ct_new);
    int m2 = LWE::Decrypt(LWEsk, ct);
    cout << m2 << "\n";

    if (message != m || message != m0 || message != m1 || message != m2){
        flag = -1;
        break;
    }
  }
  if (flag)
    cerr << "\nTests failed!\n\n";
  else
    cerr << "\nPassed all tests!\n\n";


}