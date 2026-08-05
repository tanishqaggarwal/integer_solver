#!/usr/bin/env python3
import json, re, ast, sys, glob
from propagate import NVARS
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
# precompile each equation LHS to fast eval
codes=[compile(re.sub(r'x_(\d+)', r'v[\1]', L.rsplit('=',1)[0]), '<e>','eval') for L in lines]
def count(fn):
    try: cand={int(k[2:]):v for k,v in json.load(open(fn)).items()}
    except Exception as e: return None
    v=[0]*NVARS
    for k,x in cand.items():
        if k<NVARS: v[k]=x
    ns={'v':v,'__builtins__':{}}
    f=0; fl=[]
    for i,c in enumerate(codes):
        if eval(c,ns)!=0:
            f+=1
            if len(fl)<40: fl.append(i)
    return (len(lines)-f, f, fl)
for fn in ['best/new_instance_partial_39007.json','achieved.json','achieved3.json','achieved4.json','achieved5.json']:
    r=count(fn)
    if r: print(f"{fn}: {r[0]}/{len(lines)}  ({r[1]} fail) {r[2][:20]}")
    else: print(f"{fn}: (missing)")
