

# This file was *autogenerated* from the file error_estimate.sage
from sage.all_cmdline import *   # import sage library

_sage_const_680 = Integer(680); _sage_const_2 = Integer(2); _sage_const_27 = Integer(27); _sage_const_1p4 = RealNumber('1.4'); _sage_const_512 = Integer(512); _sage_const_2048 = Integer(2048); _sage_const_7 = Integer(7); _sage_const_1p0 = RealNumber('1.0'); _sage_const_8 = Integer(8); _sage_const_1 = Integer(1); _sage_const_23 = Integer(23); _sage_const_25 = Integer(25); _sage_const_12 = Integer(12); _sage_const_20p96 = RealNumber('20.96'); _sage_const_6 = Integer(6); _sage_const_4 = Integer(4); _sage_const_3 = Integer(3)# FHEW system parameters
## LWE parameters
n = _sage_const_680 
Q = _sage_const_2 **_sage_const_27 
sigma = _sage_const_1p4 
q = _sage_const_512 
## Ring-GSW parameters
N = _sage_const_2048 
stigma = sigma
## Gadget Matrix
Bg = _sage_const_2 **_sage_const_7 
dg = ceil(_sage_const_1p0 *log(Q, Bg))
u = Q/_sage_const_8  + _sage_const_1 
## Bootstrapping Key Parameters
Br = _sage_const_23 
dr = ceil(_sage_const_1p0 *log(q, Br))
# Key Switching Key Parameters
Bks = _sage_const_25 
dks = ceil(_sage_const_1p0 *log(Q, Bks))

# standard deviation of the final error err(Refresh(c)), modelled as Gaussian distribution, the coeff 240 is debatable.
beta = sqrt((q**_sage_const_2 )/(Q**_sage_const_2 ) * (stigma**_sage_const_2  * (_sage_const_1 *Bg**_sage_const_2 )/_sage_const_12  *n*dr * (q/_sage_const_2 ) * _sage_const_2 *N*dg + sigma**_sage_const_2 *N*dks) + (n/_sage_const_2  + _sage_const_1 )/_sage_const_12 ) 
beta = sqrt(_sage_const_20p96 )
# error probability per NAND
r = (q/_sage_const_8 )/(sqrt(_sage_const_2 )*beta)
# p = int(log(1 - erf(r/sqrt(2)), 2))
x = r/sqrt(_sage_const_2 )
p = int(log(_sage_const_1 /_sage_const_6 *exp(-x*x)+_sage_const_1 /_sage_const_2 *exp(-_sage_const_4 /_sage_const_3 *x*x), _sage_const_2 )) #erfc(x) = 1-erf(x) = 1/6EXP(-x^2)+1/2EXP(-4/3*x^2), wikipedia on 'Error function'


print ("Standard deviation of NAND gate: " + str(beta))
print ("DFR per NAND gate: 2^({p})".format(p = p))

