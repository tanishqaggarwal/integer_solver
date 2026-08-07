import sys, math, random; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L, sympy
P=env.P
Np=115792089237316195423570985008687907853031073199722524052490918277602762621571
print('N\' =',Np,'prime?',sympy.isprime(Np))
# small factors
n=Np; small=[]
for q in sympy.primerange(2,3000000):
    while n%q==0: small.append(q); n//=q
    if q*q>n: break
print('small prime factors:',small)
print('cofactor:',n,'prime?',sympy.isprime(n) if n>1 else '-')
if n>1 and not sympy.isprime(n):
    print('trying pollard rho / sympy factorint (60 s budget)...',flush=True)
    try:
        f=sympy.factorint(n, limit=10**7)
        print('factorint:',f)
    except Exception as e: print('err',e)
