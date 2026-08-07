"""Find ALL OR nodes in the circuit and the true global root."""
import sys, os, collections, pickle
F='/home/user/integer_solver/solve_lab/agentT_work/mirror/F'; sys.path.insert(0,F)
from fwd import Engine, NV
from parse import node_str
E=Engine()
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
def deref(v):
    seen=set()
    while v in defrhs and defrhs[v][0]=='v' and v not in seen: seen.add(v); v=defrhs[v][1]
    return v
def as_or(v):
    v=deref(v); r=defrhs.get(v)
    if r is None or r[0]!='-': return None
    def unv(n):
        if n[0]=='v': return defrhs.get(deref(n[1]))
        return n
    Ln=unv(r[1]); Rn=unv(r[2])
    if Ln is None or Rn is None: return None
    if Ln[0]=='+' and Rn[0]=='*':
        if {node_str(Ln[1]),node_str(Ln[2])}=={node_str(Rn[1]),node_str(Rn[2])}:
            return (Ln[1],Ln[2])
    return None
ORS={}
for v in list(defrhs):
    o=as_or(v)
    if o is None: continue
    dv=deref(v)
    if dv in ORS: continue
    ch=[]
    for c in o:
        ch.append(deref(c[1]) if c[0]=='v' else None)
    ORS[dv]=tuple(ch)
print('total OR nodes in circuit:',len(ORS))
child=set()
for n,(a,b) in ORS.items(): child.add(a); child.add(b)
roots=[n for n in ORS if n not in child]
print('OR roots (not a child of any OR):',roots)
for R in roots:
    # count leaves
    def cnt(n,seen=None):
        if n not in ORS: return 1
        return cnt(ORS[n][0])+cnt(ORS[n][1])
    def dep(n):
        if n not in ORS: return 0
        return 1+max(dep(ORS[n][0]),dep(ORS[n][1]))
    lv=[]
    def col(n):
        if n not in ORS: lv.append(n); return
        col(ORS[n][0]); col(ORS[n][1])
    col(R)
    fr=[x for x in lv if x not in defrhs]
    print('  root x%d: %d leaves, depth %d, %d free / %d dead'%(R,len(lv),dep(R),len(fr),len(lv)-len(fr)))
pickle.dump(ORS,open('ors.pkl','wb'))
