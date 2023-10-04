#include <iostream>
#include <cstdlib>
#include "LWE.h"
#include "FHEW.h"
#include "FFT.h"
#include "distrib.h"
#include <time.h>

using namespace std;

LWE::SecretKeyN FHEWsk_z;
Ring_NTT FHEWsk_zNTT;
FHEW::ct_NTT FHEWctNTT;
FHEW::ct_NTT1* FHEWctNTT1;

void help(char* cmd) {
  cerr << "\nusage: " << cmd << " n\n\n" 
  << "  Generate one secret key LWEsk for LWE sample, and repeat the following test n times:\n"
  << "  prepare an LWE sample e1=LWE(m1) and an LWE sample e2=LWE(m2) under LWEsk and q/4"
  << "  perform homomorphic NAND operation to get LWE(nand(m1,m2)) under LWEsk and q/2"
  << "  then perform the bootstrapping to get LWE(nand(m1,m2)) under LWEsk and q/4"
  << "  finally validate the correctness by decrypting the LWE sample"
  << "\n If any of the tests fails, print ERROR and stop immediately.\n\n";
  exit(0);
}

double Mean(double avg_n, int x_n_plus_1, int n){
	return avg_n + (x_n_plus_1-avg_n)/(n+1); //compute avg_n_plus_1
} 

double Var(double var_n, double avg_n, double avg_n_plus_1, double x_n_plus_1, int n){
	return var_n + avg_n*avg_n - avg_n_plus_1*avg_n_plus_1 + (x_n_plus_1*x_n_plus_1 - var_n - avg_n*avg_n)/(n + 1); //compute var_n_plus_1
} 

int main(int argc, char *argv[]) {
  if (argc != 2) help(argv[0]);
  int count = atoi(argv[1]); 

  cerr << "Setting up FHEW \n";
  FHEW::Setup();
  cerr << "Generating secret key ... ";
  LWE::SecretKey LWEsk;
  LWE::KeyGen(LWEsk);
  cerr << " Done.\n";
  cerr << "Generating evaluation key ... this may take a while ... ";
  FHEW::EvalKey EK;
  FHEW::KeyGen(&EK, LWEsk);
  cerr << " Done.\n\n";
  cerr << "Testing homomorphic NAND gate " << count << " times.\n"; 

  int v1,v2;
  LWE::CipherText e1, e2, e12;
  clock_t tStart = clock();

  //statistics for DFR
  double avg = 0, avg_next = 0, var = 0;

  for (int i = 0; i < count; i++){
	v1 = rand()%2;  
  	v2 = rand()%2;
  	LWE::Encrypt(&e1, LWEsk, v1);
  	LWE::Encrypt(&e2, LWEsk, v2);

  	cerr << "Enc(" << v1 << ")";
    cerr << " NAND ";
    cerr << "Enc(" << v2 << ") = ";

    FHEW::HomNAND(&e12, EK, e1, e2);

    int v12 = LWE::Decrypt(LWEsk, e12);

    cerr << "Enc(" << v12 << ")";
    cerr << endl;

    if ((not(v1 && v2)) != v12) {
      cerr << "\n ERROR: incorrect Homomorphic Gate computation at iteration " << i+1 << "\n"; 
      exit(1); 
    }
    else { //retrieve the error value
    	int r = e12.b;
    	for (int i = 0; i < n; ++i) 
    		r -= e12.a[i] * LWEsk[i];
    	r = ((r -(q/4)*v12)%q + q)%q;
    	r = r > q/4 ? r-q : r;
    	cout << "error value: " << r << "\n";
    	//compute the mean and variance
    	avg_next = Mean(avg, r, i);
    	var = Var(var, avg, avg_next, r, i);
    	avg = avg_next;
 		cout << "mean: " << avg << "\t variance: " << var << "\n";
    }
  
  }
  cerr << "\nPassed all tests!\n\n";
  /* Do your stuff here */
  printf("Time taken per NAND gate: %.2fs\n", (double)(clock() - tStart)/(count*CLOCKS_PER_SEC));  

  return 0;
}



