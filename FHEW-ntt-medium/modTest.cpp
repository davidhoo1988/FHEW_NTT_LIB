#include <iostream>
#include <cstdlib>
#include "LWE.h"
#include "distrib.h"
#include <time.h>

using namespace std;

ZmodQ x, y, xy, xy1;


void help(char* cmd) {
  cerr << "\nusage: " << cmd << " n\n\n" 
  << "  Generate two random Z_Q elements a, b, and repeat the following test n times:\n"
  << "   - multiply a and b using either standard or proposed modular reduction \n"
  << "   - compare their results and check correctness \n"
  << "\n If any of the tests fails, print ERROR and stop immediately.\n\n";
  exit(0);
}

int main(int argc, char *argv[]) {
  if (argc != 2) help(argv[0]);
  int count = atoi(argv[1]); 
  srand(time(NULL));
  int flag = 0;
  clock_t tStart = clock();

  for (int i = 0; i < count; i++){
  	ZmodQ x = rand() % Q;
  	ZmodQ y = rand() % Q;
  	xy = LWE::ModMulQ(x, y);
    //xy1 = (x + y) % Q;
  	xy1 = uint64_t(x)*y % Q;
  	if ((xy != xy1) && (xy != xy1-Q)){
  		cout << "xy = " << xy << "\n";
   		cout << "xy1 = " << xy1 << "\n";
  		flag = -1;
  		break;
  	}
  }

  if (flag)
    cerr << "\nTests failed!\n\n";
  else
    cerr << "\nPassed all tests!\n\n";
  /* Do your stuff here */
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

}