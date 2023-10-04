# Find the bigger modulus Q s.t. q = 1 mod 2N
N = 1792
Qbit = 32
for i in range(1000000):
	if is_prime(2^Qbit-i) and (2^Qbit-i)%(2*N) == 1:
		Q = 2^Qbit-i
		print ("Q = 2^" + str(Qbit)+ "-" + str(i) + ' = ' + str(Q))
		print ("factor(Q-1) = " + str(factor(Q-1)))
		break

# Find the invertible parameter u \approx Q/2t
t = 4
F = GF(Q)
u = round(Q/(2*t))
u_inv = F(u)^(-1)
u_inv = int(u_inv)
if u_inv > Q/2:
	u_inv = u_inv - Q
print ("u = Q/2t = " + str(u))
print ("u_inverse = " + str(u_inv))

# find the primitive root of unity

gen = primitive_root(Q)
gen = mod(gen,Q)

g_list = []
prime_factor = [ i[0] for i in list(factor(Q-1)) ]

testv = [(Q-1)/i for i in prime_factor]
for gen in range(1, Q):
	flag = True
	gen = mod(gen, Q)
	for v in testv:
		if gen^v == 1:
			flag = False
			break
	if flag == True:
		# print ("primitive generator of Q is: " + str(gen))
		g_list.append(gen)
		break

root_2Nth_list = []
root_Nth_list = []
for gen in g_list:
	root_2Nth_list.append( gen^((Q-1)/(2*N)) )
	root_Nth_list.append( gen^((Q-1)/N) )

root_2Nth = root_2Nth_list[0]
root_Nth = root_Nth_list[0]

tmp = int(gen)
if tmp > Q/2:
	tmp = tmp - Q
print ( "primitive generator of Q is: " + str(tmp) )

tmp = int(root_2Nth)
if tmp > Q/2:
	tmp = tmp - Q
print ( "primitive 2N-th root of unity of Q is: " + str(tmp) )

tmp = int(root_Nth)
if tmp > Q/2:
	tmp = tmp - Q
print ( "primitive N-th root of unity of Q is: " + str(tmp) )

tmp = int(root_2Nth^(-1))
if tmp > Q/2:
	tmp = tmp - Q
print ( "inverse primitive 2N-th root of unity of Q is: " + str(tmp) )

tmp = int(root_Nth^(-1))
if tmp > Q/2:
	tmp = tmp - Q
print ( "inverse primitive N-th root of unity of Q is: " + str(tmp) )

tmp = int(mod(N, Q))
print ( "N mod Q is: " + str(tmp) )

tmp = int( (mod(N, Q))^(-1) )
if tmp > Q/2:
	tmp = tmp - Q
print ( "inverse N mod Q is: " + str(tmp) )


# find omega's psi's for NTT computation
omega_list = []
psi_list = []
omega_inv_list = []
psi_inv_list = []
for i in range(N):
	omega_list.append(root_Nth^i)
	psi_list.append(root_2Nth^i)
	omega_inv_list.append(root_Nth^(-i))
	psi_inv_list.append(root_2Nth^(-i))

# find vgprime = [u, u<<6, u<<12, u<<18, u<<24, u<<30]
vgprime = []
Bg = 2^6
log_Bg = log(Bg, 2)
for i in range(ceil(Qbit/log_Bg)):
	tmp = int(mod(u*Bg^(i),Q))
	if tmp > Q/2:
		tmp = tmp - Q
	vgprime.append(tmp)
print (vgprime)

# estimate the security by using uSVP, decoding, and dual attacks.
n = 505
q = 256
sigma = 1.4
load("../../../malb-lwe-estimator-fb7deba98e59/estimator.py")

# estimate LWE
print ("LWE security esitmation as follows:")
n = n; q = q; stddev = sigma; alpha = alphaf(sigmaf(stddev), q)
_ = estimate_lwe(n, alpha, q, secret_distribution=((-1,1),n/2), reduction_cost_model=BKZ.sieve)

# estimate RLWE
print ("RLWE security esitmation as follows:")
n = N; q = Q; stddev = sigma; alpha = alphaf(sigmaf(stddev), q)
_ = estimate_lwe(n, alpha, q, reduction_cost_model=BKZ.sieve)