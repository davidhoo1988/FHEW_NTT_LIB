#include <iostream>
#include <cstdlib>
#include "LWE.h"
#include "FHEW.h"
#include "FFT.h"
#include "distrib.h"
#include <time.h>

using namespace std;

Ring_ModQ a, b, ab, ab1, absum, ab1sum, ab2sum, tmp;
Ring_FFT aFFT, bFFT, abFFT, absumFFT;
Ring_NTT aNTT, bNTT, abNTT, absumNTT;

void help(char* cmd) {
  cerr << "\nusage: " << cmd << " n\n\n" 
  << "  Generate two random ring elements a, b (one of them have coeff smaller than B_g), and repeat the following test n times:\n"
  << "   - multiply a and b using either FFT or anti-cyclic conovlution\n"
  << "   - compare their results and check correctness \n"
  << "\n If any of the tests fails, print ERROR and stop immediately.\n\n";
  exit(0);
}

int main(int argc, char *argv[]) {
  if (argc != 2) help(argv[0]);
  int count = atoi(argv[1]); 
  srand(time(NULL));
  int flag = 0;

  FFTsetup();
  NTTsetup();
  clock_t tStart = clock();

  for (int test = 0; test < count; test++){
    for (int i = 0; i < K2; i++){
        //cerr << "Setting up operand a \n";
        for (int i = 0; i < N; i++)
        	a[i] = LWE::ModMulQ(rand(),1); //rand()%(1<<16) 
        // FFTforward(aFFT, a);
        // FFTbackward(tmp, aFFT);
        // for (int i = 0; i < N; i++)
      	 //  if (a[i] != tmp[i]){
      	 //  	cout << "FFT(a) incorrect!\n";
      	 //  	flag = -1;
      	 //  	return -1;
      	 //  }
        NTTforward(aNTT, a);
        NTTbackward(tmp, aNTT);
        for (int i = 0; i < N; i++)
          if (a[i] != tmp[i]){
            cout << "NTT(a) incorrect!\n";
            flag = -1;
            return -1;
          }

        //cerr << "Setting up operand b \n";
        for (int i = 0; i < N; i++)
        	b[i] = LWE::ModMulQ(rand(),1); //(1<<5);
        // FFTforward(bFFT, b);
        // FFTbackward(tmp, bFFT);
        // for (int i = 0; i < N; i++)
      	 //  if (b[i] != tmp[i]){
      	 //  	cout << "FFT(b) incorrect!\n";
      	 //  	flag = -1;
      	 //  	return -1;
      	 //  } 
        NTTforward(bNTT, b);
        NTTbackward(tmp, bNTT);
        for (int i = 0; i < N; i++)
          if (b[i] != tmp[i]){
            cout << "NTT(b) incorrect!\n";
            flag = -1;
            return -1;
          }
        //cerr << "Computing FFT\n";    
        AntiCyclicConv(ab, a,b);
        if (i == 0)
          for (int j = 0; j < N; j++)
            absum[j] = ab[j];
        else 
          for (int j = 0; j < N; j++)
            absum[j] = LWE::ModAddQ(absum[j], ab[j]);  

        // for (int j = 0; j < N2; j++)
        //   abFFT[j] = aFFT[j]*bFFT[j];
        // if (i == 0)
        //   for (int j = 0; j < N2; j++)
        //     absumFFT[j] = abFFT[j];
        // else 
        //   for (int j = 0; j < N2; j++)
        //     absumFFT[j] += abFFT[j];  

        for (int j = 0; j < N; j++)
          abNTT[j] = LWE::ModMulQ(aNTT[j], bNTT[j]);
        if (i == 0)
          for (int j = 0; j < N; j++)
            absumNTT[j] = abNTT[j];
        else 
          for (int j = 0; j < N; j++)
            absumNTT[j] = LWE::ModAddQ(absumNTT[j], abNTT[j]);  
    }

    // FFTbackward(ab1sum, absumFFT);
    NTTbackward(ab2sum, absumNTT);
    for (int k = 0; k < N; k++)
    if (absum[k] != ab2sum[k]){ //absum[k] != ab2sum[k]
      cout << "test = " << k << "\n";
      cout << "absum[k] = " << absum[k] << "\n";
      cout << "ab2sum[k] = " << ab2sum[k] << "\n";
      flag = -1;
    }
    if (flag)
        break;
  }

  NTTdestroy();
  if (flag)
    cerr << "\nTests failed!\n\n";
  else
    cerr << "\nPassed all tests!\n\n";
  /* Do your stuff here */
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  return 0;
}