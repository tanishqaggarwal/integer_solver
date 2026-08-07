"""Is the map (message bits) -> (A,B) AFFINE over F_p?  Test additivity on pairs."""
import os, sys, itertools, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import g29_frame as F
import gsym2 as G
from gsym2 import L, ad, P
c11,c12,c21,c22=8646263,1073965,10159099,6926539
det=(c11*c22-c12*c21)%P; inv=pow(det,-1,P)
def AB(FL):
    v=list(F.v0)
    for b in FL: v[b]=1-v[b]
    ad.fwd(v,rounds=8)
    r=F.analyse(v)
    if r['incchecks'] or r['nzc']: return None, r
    d={a:(g%P if isinstance(g,int) else None) for a,g in r['res']}
    if 19297 not in d or 19299 not in d: 
        # residual may be zero for those atoms
        d.setdefault(19297,0); d.setdefault(19299,0)
    if d[19297] is None or d[19299] is None: return 'POLY', r
    a1,a2=d[19297],d[19299]
    A=(c22*a1-c12*a2)*inv%P
    B=(-c21*a1+c11*a2)*inv%P
    return (A,B), r
base,_=AB([])
print('base (A,B) =',base)
random.seed(5)
BITS=[47,91,112,438,490,542,853,1203,1357,1413]
vals={}
for b in BITS:
    ab,r=AB([b])
    vals[b]=ab
    print('bit x%-6d (A,B)=%s  nres=%d'%(b,str(ab)[:80],len(r['res'])))
print('--- pair additivity test ---')
for i,j in itertools.combinations(BITS[:5],2):
    ab,r=AB([i,j])
    if not isinstance(ab,tuple) or not isinstance(vals[i],tuple) or not isinstance(vals[j],tuple):
        print('  (%d,%d): non-numeric'%(i,j)); continue
    pred=((vals[i][0]+vals[j][0]-base[0])%P,(vals[i][1]+vals[j][1]-base[1])%P)
    print('  (%d,%d): measured=%s predicted_affine=%s  %s'%(i,j,str(ab)[:60],str(pred)[:60],'AFFINE' if ab==pred else 'NOT affine'))
