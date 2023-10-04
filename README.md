# FHEW-ntt library

## basic information on FHEW-ntt library
FHEW-ntt is an NTT-based version FHEW released here: https://github.com/lducas/FHEW. 
Read the classic paper "FHEW: bootstrapping homomorphic encryption in less than a second"  to undertand more on FHEW.
The NTT implementation is adpated from Nayuki-FFT-Portable, https://www.nayuki.io/page/free-small-fft-in-multiple-languages.
Also ensure fftw library (see https://www.fftw.org/download.html) is installed before you proceed to install FHEW-ntt library.


## subdirectories in FHEW-ntt
FHEW-ntt-128 is an NTT-based FHEW with 128-bit security level, read ./FHEW-ntt-128/README.md to learn how to compile and install it.

FHEW-ntt-192 is an NTT-based FHEW with 192-bit security level, read ./FHEW-ntt-192/README.md to learn how to compile and install it. 

FHEW-ntt-256 is an NTT-based FHEW with 256-bit security level, read ./FHEW-ntt-256/README.md to learn how to compile and install it.


## Other remarks
FHEW-ntt-medium and  FHEW-ntt-128AP are identical to FHEW-ntt-128, and is kept for historical reasons.  'medium' suggests the security level is medium if compared with the 80-bit FHEW-ntt. 'AP' suggests the external product is asymmetric, and is performed between a RLWE ciphertext and a RGSW ciphertext for better computational efficiency.
