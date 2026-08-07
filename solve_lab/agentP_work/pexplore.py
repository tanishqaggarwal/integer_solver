#!/usr/bin/env python3
import pickle, json, sys, re
from collections import defaultdict, Counter
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb'))
rows,AP=D['rows'],D['AP']
v=[0]*38748
for k,val in json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')).items():
    v[int(k[2:])]=int(val)

var2at=defaultdict(list)
for i,ap in enumerate(AP):
    s=set()
    for m in ap: s.update(m)
    for x in s: var2at[x].append(i)

def fmt(ap):
    ts=[]
    for m,c in sorted(ap.items()):
        cs=('' if c==1 else '-' if c==-1 else str(c)+'*') if m else str(c)
        if len(str(abs(c)))>12: cs='K%d*'%(abs(c)%100000) if m else 'K'
        ts.append(cs+'*'.join('x%d'%i for i in m) if m else cs)
    return ' + '.join(ts).replace('+ -','- ')

def ev(ap):
    s=0
    for m,c in ap.items():
        t=c
        for i in m: t*=v[i]
        s+=t
    return s

if __name__=='__main__':
    vals=[ev(a) for a in AP]
    nz=[i for i,x in enumerate(vals) if x]
    print("nonzero atoms:",len(nz))
    bad=[]
    for ei,r in enumerate(rows):
        s=sum(c*vals[a] for c,a in r['row'])
        if r['scal']*s**r['pw']!=0: bad.append(ei)
    print("failing eqs:",bad)
    pickle.dump({'vals':vals},open(W+'atomvals4.pkl','wb'))
    for a in (14823,17423,15534,37280,37985):
        print(f"--- x{a} (val bits {v[a].bit_length()}) in {len(var2at[a])} atoms:")
        for i in var2at[a][:12]:
            print("    ",fmt(AP[i]), " =0?", vals[i]==0)
