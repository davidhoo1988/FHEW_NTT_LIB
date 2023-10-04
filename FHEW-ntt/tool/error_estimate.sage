# FHEW system parameters
## LWE parameters
n = 500
Q = 2^32
sigma = 2^17
q = 2^9
## Ring-GSW parameters
N = 2^10
stigma = 3.2
## Gadget Matrix
Bg = 2^6
dg = 6
u = Q/8 + 1
## Bootstrapping Key Parameters
Br = 23
dr = 2
# Key Switching Key Parameters
Bks = 24
dks = 7

# standard deviation of the final error err(Refresh(c)), modelled as Gaussian distribution
beta = sqrt((q^2)/(Q^2) * (stigma^2 * (2400*Bg^2)/12 *n*dr * (q/2) * 2*N*dg + sigma^2*N*dks) + (n/2 + 1)/12)
# beta = 6.94
# error probability per NAND
r = (q/8)/(sqrt(2)*beta)
p = int(log(1 - erf(r/sqrt(2)), 2))
# p = 1 - erf(r/sqrt(2))

print (beta)
print (p)