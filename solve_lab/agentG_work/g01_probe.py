import os, sys, json
LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff'))
sys.path.insert(0, os.path.join(LAB,'s10'))
import lib as L
P = 2**256-2**32-977
src = sys.argv[1] if len(sys.argv)>1 else os.path.join(LAB,'best','new_instance_partial_39026.json')
v = L.load(src)
av = L.all_atom_values(v)
fail = L.failing_eqs(av)
print('state', src, 'score', L.NEQ-len(fail), 'failing', fail[:20])
coords = {'x1':12186,'y1':16742,'x2':14853,'y2':24908,'x3':22162,'y3':30213}
vals = {k:(v[i]%P) for k,i in coords.items()}
for k,i in coords.items():
    print(f'  {k}=x{i}: raw bits {v[i].bit_length()}  r={v[i]%P}')
def oncurve(x,y): return (y*y-x*x*x-7)%P
for pt in [('P1','x1','y1'),('P2','x2','y2'),('P3','x3','y3')]:
    n,xa,ya=pt
    print(f'  {n} on curve y^2=x^3+7 ? residual={oncurve(vals[xa],vals[ya])}')
# addition check
x1,y1,x2,y2,x3,y3=[vals[k] for k in ['x1','y1','x2','y2','x3','y3']]
if (x2-x1)%P:
    lam=(y2-y1)*pow(x2-x1,-1,P)%P
    ax=(lam*lam-x1-x2)%P; ay=(lam*(x1-(lam*lam-x1-x2))-y1)%P
    print('  P1+P2 =', ax, ay)
    print('  matches (x3,y3)?', ax==x3%P, ay==y3%P)
