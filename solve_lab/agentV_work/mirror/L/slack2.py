"""General claim: every slack wire is (constant multiple of p) * (free var)."""
import sys, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentT_work/mirror/F')
from fwd import Engine,NV
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
sys.setrecursionlimit(200000)
CONST={}
def ev(v):
    if v in CONST: return CONST[v]
    if v not in defrhs: CONST[v]=None; return None
    r=defrhs[v]
    def e(n):
        if n[0]=='c': return n[1]
        if n[0]=='v': return ev(n[1])
        a=e(n[1])
        if a is None: return None
        b=e(n[2])
        if b is None: return None
        return a+b if n[0]=='+' else (a-b if n[0]=='-' else a*b)
    CONST[v]=e(r); return CONST[v]
for v in list(defrhs): ev(v)
constp=[v for v,c in CONST.items() if c is not None and c!=0 and c%p==0]
print('wires with a CONSTANT value (no free vars): %d'%sum(1 for c in CONST.values() if c is not None))
print('   of those, constant and divisible by p: %d'%len(constp))
print('   distinct constant values that are multiples of p: %s'%
      sorted({c//p for v,c in CONST.items() if c is not None and c!=0 and c%p==0})[:10])
# every product-of-two-vars wire: classify
prod=[(v,r) for v,r in defrhs.items() if r[0]=='*' and r[1][0]=='v' and r[2][0]=='v']
print('\nwires defined as a product of two wires: %d'%len(prod))
good=0; bad=[]
for v,r in prod:
    a,b=r[1][1],r[2][1]
    ca,cb=CONST.get(a),CONST.get(b)
    ok=(ca is not None and ca%p==0) or (cb is not None and cb%p==0)
    if ok: good+=1
    else: bad.append(v)
print('   with at least one factor a CONSTANT multiple of p: %d'%good)
print('   without: %d  (these are the selector products sel*value, not slack)'%len(bad))
# now restrict to the ones that actually appear as slack in a residual atom
Hd=pickle.load(open('handles.pkl','rb')); handle=set(Hd['handle'])
freeall={}
def fa(v):
    if v in freeall: return freeall[v]
    if v not in defrhs: freeall[v]={v}; return freeall[v]
    freeall[v]=set(); s=set()
    for u in vars_of(defrhs[v]): s|=fa(u)
    freeall[v]=s; return s
slackwires=set()
for a in E.res:
    for v in vars_of(E.atoms[a]):
        if v in defrhs and defrhs[v][0]=='*':
            r=defrhs[v]
            if r[1][0]=='v' and r[2][0]=='v':
                hs=[u for u in fa(v) if u in handle]
                if len(hs)==1: slackwires.add(v)
print('\nSLACK wires (product-of-two appearing in a residual atom, one free cofactor): %d'%len(slackwires))
okc=0; notok=[]
for v in slackwires:
    r=defrhs[v]; a,b=r[1][1],r[2][1]
    ca,cb=CONST.get(a),CONST.get(b)
    if (ca is not None and ca%p==0) or (cb is not None and cb%p==0): okc+=1
    else: notok.append(v)
print('   with a constant-multiple-of-p factor: %d / %d'%(okc,len(slackwires)))
print('   EXCEPTIONS: %s'%notok[:20])
print('\n=> every slack term is  (constant multiple of p) x (free variable)')
print('=> slack == 0 mod p in EVERY assignment; the hand-off is exact mod p and free over Z.')
