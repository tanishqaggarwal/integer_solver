#!/usr/bin/env python3
"""Agent P: are the 383 repeating 43-atom blocks structurally identical?"""
import pickle,sys,json
from collections import Counter,defaultdict
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']
S=pickle.load(open(W+'slp.pkl','rb')); topo=S['topo']; outof=S['outof']
Qv=24453
qpos=[i for i in range(len(topo)) if any(Qv in m for m in AP[topo[i]]) and len(AP[topo[i]])==3]
print("Q gates:",len(qpos),"spacing set:",set(b-a for a,b in zip(qpos,qpos[1:])))
LO,HI=-6,37   # block window relative to Q gate

def blocksig(q, keepcoef=True):
    ren={}; out=[]
    for i in range(q+LO,q+HI+1):
        if i<0 or i>=len(topo): return None
        ap=AP[topo[i]]
        terms=[]
        for m,c in sorted(ap.items(),key=lambda z:(len(z[0]),z[0])):
            mm=[]
            for x in m:
                if x==Qv: mm.append('Q')
                else:
                    if x not in ren: ren[x]=len(ren)
                    mm.append(ren[x])
            cc = c if keepcoef else ('C' if abs(c)>2 else c)
            if abs(c)>10**20: cc='BIG'
            terms.append((tuple(sorted(mm,key=str)),cc))
        out.append(tuple(sorted(terms,key=str)))
    return tuple(out), ren

sigs_c=Counter(); sigs_n=Counter(); byn=defaultdict(list)
for q in qpos:
    r=blocksig(q,True); sigs_c[r[0]]+=1
    r2=blocksig(q,False); sigs_n[r2[0]]+=1; byn[r2[0]].append(q)
print("distinct block signatures WITH coefficients:",len(sigs_c))
print("distinct block signatures IGNORING small coefficients:",len(sigs_n))
for s,c in sigs_n.most_common(6): print("   count",c,"first q:",byn[s][:4])
pickle.dump({'qpos':qpos},open(W+'qpos.pkl','wb'))
