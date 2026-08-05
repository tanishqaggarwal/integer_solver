import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
import sy_build as B
p=H.p
B.regime11()
base=H.val[:]
def loads():
    return H.val[31731], H.val[9106], H.val[2239]
# linear model f(x8731,x9118) = f0 + c8*x8731 + c9*x9118
def setk(a,b):
    H.val[8731]=a; H.val[9118]=b; H.forward()
setk(base[8731],base[9118]); f31_0,f91_0,f22_0 = loads()
# derivative wrt x8731
setk(base[8731]+1,base[9118]); f31_a,f91_a,f22_a=loads()
# derivative wrt x9118
setk(base[8731],base[9118]+1); f31_b,f91_b,f22_b=loads()
setk(base[8731],base[9118])
c31_8=f31_a-f31_0; c31_9=f31_b-f31_0
c91_8=f91_a-f91_0; c91_9=f91_b-f91_0
c22_8=f22_a-f22_0; c22_9=f22_b-f22_0
# absolute constants (value at x8731=x9118=0):
A31_0=f31_0 - c31_8*base[8731] - c31_9*base[9118]
A91_0=f91_0 - c91_8*base[8731] - c91_9*base[9118]
A22_0=f22_0 - c22_8*base[8731] - c22_9*base[9118]
print('x_31731 = %d + (%d)*x8731 + (%d)*x9118'%(A31_0,c31_8,c31_9))
# verify linearity at base
def chk(a,b):
    setk(a,b); v=loads(); setk(base[8731],base[9118])
    pred=(A31_0+c31_8*a+c31_9*b, A91_0+c91_8*a+c91_9*b, A22_0+c22_8*a+c22_9*b)
    return v==pred
print('linear model exact at base+random?', chk(base[8731]+12345,base[9118]-999))
# Solve x_31731=0: A31_0 + c31_8*x8 + c31_9*x9 = 0  (Diophantine)
import math
g=math.gcd(c31_8,c31_9)
print('gcd(c31_8,c31_9)=',g,' A31_0 % g ==0?', A31_0%g==0)
# general solution: x8 = x8p + (c31_9//g)*t ; x9 = x9p - (c31_8//g)*t
# particular via extended euclid
def egcd(a,b):
    if b==0: return (a,1,0)
    gg,x,y=egcd(b,a%b); return (gg,y,x-(a//b)*y)
gg,s,t0=egcd(c31_8,c31_9)   # c31_8*s + c31_9*t0 = gg
assert gg==g or gg==-g
# want c31_8*x8 + c31_9*x9 = -A31_0
k=-A31_0//gg
x8p=s*k; x9p=t0*k
# verify
assert A31_0 + c31_8*x8p + c31_9*x9p==0
dx8=c31_9//gg; dx9=-c31_8//gg
print('particular x8p,x9p ok; direction (dx8,dx9)=(%d,%d)'%(dx8,dx9))
# Now x_9106(t) and x_2239(t) along the family:
# x9106 = A91_0 + c91_8*(x8p+dx8*t) + c91_9*(x9p+dx9*t) = base91 + slope91*t
base91=A91_0+c91_8*x8p+c91_9*x9p; slope91=c91_8*dx8+c91_9*dx9
base22=A22_0+c22_8*x8p+c22_9*x9p; slope22=c22_8*dx8+c22_9*dx9
print('x9106(t)=%d + %d*t'%(base91,slope91))
print('x2239(t)=%d + %d*t'%(base22,slope22))
# Condition (3): p | x2239(t): base22 + slope22*t ≡0 mod p -> t ≡ -base22/slope22 mod p
# Condition (2): 13523997*p | x9106(t)
import pickle
pickle.dump(dict(x8p=x8p,x9p=x9p,dx8=dx8,dx9=dx9,base91=base91,slope91=slope91,base22=base22,slope22=slope22,
                 base8731=base[8731],base9118=base[9118]),open('sy_dio.pkl','wb'))
# solve (3) mod p
def inv(x,m): return pow(x%m,-1,m)
M2=p; M9=13523997*p
t3 = (-base22 * inv(slope22,M2)) % M2 if slope22%p!=0 else None
print('cond3 (p|x2239): slope22%p=',slope22%p,' t ≡',t3,'mod p' if t3 is not None else 'slope22≡0!')
t2 = (-base91 * inv(slope91,M9)) % M9 if math.gcd(slope91,M9)==1 else ('gcd',math.gcd(slope91,M9))
print('cond2 (13523997p|x9106): gcd(slope91,M9)=',math.gcd(slope91,M9))
