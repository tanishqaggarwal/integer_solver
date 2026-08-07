#!/usr/bin/env python3
import pickle,sys,json
from collections import defaultdict,Counter
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb'))
rows,AP=D['rows'],D['AP']
S=pickle.load(open(W+'slp.pkl','rb'))
topo=S['topo']; outof=S['outof']
gold=[0]*38748
for k,val in json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')).items():
    gold[int(k[2:])]=int(val)

def K(c):
    a=abs(c)
    if a>10**20: return ('-' if c<0 else '')+'<%d-bit>'%a.bit_length()
    return str(c)
def fmt(ap,o=None):
    ts=[]
    for m,c in sorted(ap.items(),key=lambda z:(len(z[0]),z[0])):
        if not m: ts.append(K(c)); continue
        mono='*'.join('x%d'%i for i in m)
        if c==1: ts.append(mono)
        elif c==-1: ts.append('-'+mono)
        else: ts.append(K(c)+'*'+mono)
    s=' + '.join(ts).replace('+ -','- ')
    return (f"[def x{o}] " if o is not None and o>=0 else "[CONSTR ] ")+s+" = 0"

def show(lo,hi):
    for i in range(lo,hi):
        a=topo[i]
        print(f"{i:6d} {fmt(AP[a],outof[a])}")

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='find':
        # locate leaf-load atoms in topo order
        big=[]
        for i,a in enumerate(topo):
            if any(abs(c)>10**20 for c in AP[a].values()): big.append(i)
        print("leaf-load atom positions:",len(big))
        print(big[:40])
        import numpy as np
        d=[big[j+1]-big[j] for j in range(len(big)-1)]
        print("gap hist:",sorted(Counter(d).items())[:20])
    else:
        show(int(sys.argv[1]),int(sys.argv[2]))
