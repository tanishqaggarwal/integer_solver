#!/usr/bin/env python3
import pickle,sys,json
from collections import defaultdict,Counter
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb'))
AP=D['AP']
S=pickle.load(open(W+'slp.pkl','rb'))
topo=S['topo']; outof=S['outof']
g=[0]*38748
for k,v in json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')).items(): g[int(k[2:])]=int(v)
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
Q=97553848499418123410591666447050222001188385549510401465815187079080512838891
def tag(x):
    v=g[x]
    if v==P: return 'P'
    if v==Q: return 'Q'
    if v==0: return ''
    if v==1: return '1'
    return '#%d'%v.bit_length()
def K(c):
    a=abs(c)
    if a>10**20: return ('-' if c<0 else '')+'<%dbit>'%a.bit_length()
    return str(c)
def fmt(ap,o=None):
    ts=[]
    for m,c in sorted(ap.items(),key=lambda z:(len(z[0]),z[0])):
        if not m: ts.append(K(c)); continue
        mono='*'.join('x%d%s'%(i,('{'+tag(i)+'}' if tag(i) else '')) for i in m)
        ts.append(mono if c==1 else '-'+mono if c==-1 else K(c)+'*'+mono)
    s=' + '.join(ts).replace('+ -','- ')
    return (f"x{o}:= " if o is not None and o>=0 else "CONSTR: ")+s+" = 0"
def main():
  for i in range(int(sys.argv[1]),int(sys.argv[2])):
    a=topo[i]
    print(f"{i:6d} {fmt(AP[a],outof[a])}")

if __name__=="__main__": main()
