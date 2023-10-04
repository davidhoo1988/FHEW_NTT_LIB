# This script construts the cumulative density table for Gaussian sampling
def rho(x, sigma):
    return exp(-(x/sigma)^2/2)

sigma = 1.4 # set up the stdev of the Gaussian distribution
offset = 11 # set up the domain of the Gaussain dist. to [-offset,...,0,...,offset]
S = 0
for i in range(-offset,offset+1):
	S += rho(i, sigma)

cdt = [0]*(2*offset+1);
cdt[0] = rho(-offset,sigma)/S
for i in range(-offset+1,offset+1):
	cdt[i+offset] = cdt[i+offset-1] + rho(i,sigma)/S 

# print out every entry in the CDT table
for i in range(-offset,offset+1):
	print (float(cdt[i+offset]))