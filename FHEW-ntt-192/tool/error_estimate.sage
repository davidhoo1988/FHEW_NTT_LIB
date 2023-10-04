# FHEW system parameters
## LWE parameters
n = 512
Q = 2^27
sigma = 1.4
q = 2^9
## Ring-GSW parameters
N = 2^10
stigma = sigma
## Gadget Matrix
Bg = 2^7
dg = ceil(1.0*log(Q, Bg))
u = Q/8 + 1
## Bootstrapping Key Parameters
Br = 23
dr = ceil(1.0*log(q, Br))
# Key Switching Key Parameters
Bks = 25
dks = ceil(1.0*log(Q, Bks))

# standard deviation of the final error err(Refresh(c)), modelled as Gaussian distribution, the coeff 240 is debatable.
beta = sqrt((q^2)/(Q^2) * (stigma^2 * (1*Bg^2)/12 *n*dr * (q/2) * 2*N*dg + sigma^2*N*dks) + (n/2 + 1)/12) 
beta = sqrt(22.3)
# error probability per NAND
r = (q/8)/(sqrt(2)*beta)
# p = int(log(1 - erf(r/sqrt(2)), 2))
# p = 1 - erf(r/sqrt(2))
x = r/sqrt(2)
p = int(log(1/6*exp(-x*x)+1/2*exp(-4/3*x*x), 2)) #erfc(x) = 1-erf(x) = 1/6EXP(-x^2)+1/2EXP(-4/3*x^2), wikipedia on 'Error function'

print ("Standard deviation of NAND gate: " + str(beta))
print ("DFR per NAND gate: 2^({p})".format(p = p))