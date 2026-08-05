import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
import sy_build as B
p=H.p
B.regime11()
base=H.val[:]
# linear model of x_9106, x_2239 (mod p) in x_8731,x_9118
def loads():
    return H.val[9106]%p, H.val[2239]%p
def setk(a,b): H.val[8731]=a; H.val[9118]=b; H.forward()
setk(0,0); f91_0,f22_0=loads()
setk(1,0); f91_a,f22_a=loads()
setk(0,1); f91_b,f22_b=loads()
setk(base[8731],base[9118])
c91_8=(f91_a-f91_0)%p; c91_9=(f91_b-f91_0)%p
c22_8=(f22_a-f22_0)%p; c22_9=(f22_b-f22_0)%p
# solve: f91_0 + c91_8*x8 + c91_9*x9 ≡0 ; f22_0 + c22_8*x8 + c22_9*x9 ≡0
det=(c91_8*c22_9 - c91_9*c22_8)%p
def inv(x): return pow(x%p,-1,p)
r1=(-f91_0)%p; r2=(-f22_0)%p
x8=((r1*c22_9 - r2*c91_9)*inv(det))%p
x9=((c91_8*r2 - c22_8*r1)*inv(det))%p
print('det!=0?',det!=0)
print('x_8731 ≡',x8,'mod p')
print('x_9118 ≡',x9,'mod p')
# TEST with these small residues; keep G1,G2 via x_4432=x_8731, x_7068=x_9118
H.val[8731]=x8; H.val[9118]=x9
H.val[9413]=0; H.val[17325]=0
H.forward()
H.val[4432]=H.val[19964]   # G2: x_4432 = x_19964(=x_8731)
H.val[7068]=H.val[2099]    # G1: x_7068 = x_2099(=x_9118) (since x_642=0)
H.forward()
# check loads mod p
print('x_9106%p=',H.val[9106]%p,'x_2239%p=',H.val[2239]%p,'x_31731%p=',H.val[31731]%p)
# set handles
H.val[950]=H.val[9106]//(13523997*p)
H.val[6947]=(6122989*H.val[2239])//p
H.val[33168]=-(H.val[31731]//p)
H.forward()
F=H.fails()
G1=7376877*H.val[642]+H.val[2099]-H.val[7068]; G2=H.val[4432]-H.val[19964]-H.val[28730]
print('G1=',G1,'G2=',G2)
print('FAILS=',len(F)); print(sorted(F))
