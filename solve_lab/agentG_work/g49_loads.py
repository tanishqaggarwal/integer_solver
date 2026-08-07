"""Structural extraction of the load-pin table: for every boolean free input b, the
constants HUGE it pins, and whether HUGE mod p is a valid x-coordinate of the curve."""
import os, sys, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gpt, gGclose
import gsym2 as G
from gsym2 import L, ad, P
K=gpt.K; inv3=pow(3,-1,P); IU2=gpt.IU2
v=L.load('/home/user/integer_solver/solve_lab/s10/AG_39013.json'); ad.fwd(v,rounds=6)
FREE=[u for u in range(L.NVARS) if u not in L.definer]
BOOL=[u for u in FREE if gGclose.isbool(u)]
def isxcoord(r):
    X=((r+K*inv3)%P)*IU2%P
    t=(pow(X,3,P)+7)%P
    return t==0 or pow(t,(P-1)//2,P)==1
tab=collections.defaultdict(list)
nload=0
for b in BOOL:
    for a in L.var_atoms[b]:
        pl=L.polys[a]
        lin=pl.get((b,))
        if lin is None or abs(lin)<10**30: continue
        for m,c in pl.items():
            if len(m)==2 and b in m:
                xB=m[0] if m[1]==b else m[1]
                if xB==b: continue
                HUGE=(-lin*pow(c,-1,P))%P
                tab[b].append((a,xB,HUGE))
                nload+=1
print('boolean free inputs with a huge load pin: %d  (total pins %d)'%(len(tab),nload))
ok=0; tot=0
xs={}
for b,ls in tab.items():
    for a,xB,H in ls:
        tot+=1
        if isxcoord(H): ok+=1; xs.setdefault(b,[]).append((xB,H))
print('load constants whose residue is a VALID x-coordinate of the curve: %d of %d (random expectation ~%.0f)'%(ok,tot,tot/2))
print('bits with at least one valid-x load: %d'%len(xs))
pickle.dump({'tab':dict(tab),'validx':xs},open('/home/user/integer_solver/solve_lab/agentG_work/loads.pkl','wb'))
for b in sorted(xs)[:10]:
    print('  x%-6d -> %s'%(b,[(xB,str(H)[:24]+'..') for xB,H in xs[b]][:3]))
