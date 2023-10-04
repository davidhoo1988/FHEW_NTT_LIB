#ifndef FFT_H
#define FFT_H

#include "params.h"

void FFTsetup();
void FFTforward(Ring_FFT res, const Ring_ModQ val);
void FFTbackward(Ring_ModQ res, const Ring_FFT val);
void AntiCyclicConv(Ring_ModQ val, const Ring_ModQ op1, const Ring_ModQ op2);

void NTTsetup();
void NTTforward(Ring_NTT res, const Ring_ModQ val);
void NTTbackward(Ring_ModQ res, const Ring_NTT val);
void NTTdestroy();
#endif
