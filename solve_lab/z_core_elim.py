import os,sys,json
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
# CORE: S = x_33469*x_29322^2 - x_3558^2, T = x_27713*x_29322 - x_3558*x_1326
# Eliminate x_29322, x_3558 (both 0 at agentA). Resultant/elimination:
# If x_29322!=0: r=x_3558/x_29322 => x_33469=r^2, x_27713=r*x_1326 => x_27713^2 = x_33469*x_1326^2
x33469=H.val[33469]%p; x27713=H.val[27713]%p; x1326=H.val[1326]%p
x29322=H.val[29322]%p; x3558=H.val[3558]%p
print(f"agentA core vars mod p: x_29322={x29322}, x_3558={x3558} (degenerate: both {'0' if x29322==0==x3558 else 'nonzero'})")
elim = (x27713*x27713 - x33469*x1326*x1326)%p
print(f"CORE elimination ideal relation  x_27713^2 - x_33469*x_1326^2 mod p = {elim}")
print(f"  => {'agentA LIES on non-degenerate core variety (QR consistent)' if elim==0 else 'agentA does NOT lie on non-degenerate variety (strictly degenerate)'}")
# also S,T at agentA
S=(x33469*x29322*x29322 - x3558*x3558)%p
T=(x27713*x29322 - x3558*x1326)%p
print(f"S mod p={S}, T mod p={T} (both 0 => core satisfied)")
# QR of x_33469
leg=pow(x33469,(p-1)//2,p)
print(f"x_33469 Legendre = {'QR' if leg==1 else 'NQR' if leg==p-1 else '0'}")
