#!/usr/bin/env python3
"""Agent P: atom evaluation under the deliverable, + definition-DAG extraction."""
import pickle, json, sys
from collections import defaultdict, Counter
sys.set_int_max_str_digits(10_000_000)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D = pickle.load(open(W+'model2.pkl','rb'))
rows, AP = D['eq_rows'], D['atom_polys']

v=[0]*38748
for k,val in json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')).items():
    v[int(k[2:])]=int(val)

def ev(ap):
    s=0
    for m,c in ap.items():
        t=c
        for i in m: t*=v[i]
        s+=t
    return s

vals=[ev(a) for a in AP]
nz=[i for i,x in enumerate(vals) if x]
print("nonzero atoms under deliverable:", len(nz), "of", len(AP))
for i in nz[:20]:
    print("  atom",i,"val bits",vals[i].bit_length(), "poly", {k:(c if abs(c)<10**8 else str(c)[:12]+'..') for k,c in AP[i].items()})

# which equations fail
badeq=[]
for ei,r in enumerate(rows):
    s=sum(c*vals[a] for c,a in r['row'])
    if r['scal']*(s**r['pw'])!=0: badeq.append(ei)
print("failing equations from my model:", badeq)
pickle.dump({'vals':vals,'nz':nz}, open(W+'atomvals.pkl','wb'))
