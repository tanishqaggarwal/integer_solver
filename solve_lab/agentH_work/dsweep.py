"""Detach-set sweep: ALL singletons, ALL pairs, PRUNED triples.  Criterion: deficit<4 or rank>7."""
import json, time, itertools, sys
from detach import make, price, POOL, WIT
def key(r): return (r['deficit']-r['rank'], -r['score'])
res=[]; t0=time.time()
cal=price(make(WIT),'CAL'); print('CAL score=%d deficit=%d rank=%d'%(cal['score'],cal['deficit'],cal['rank']),flush=True)
CALGAP=cal['deficit']-cal['rank']
wins=[]
def rec(D,r):
    r=dict(r); r['D']=list(D); res.append(r)
    if r['deficit']<4 or r['rank']>7 or r['score']>39026:
        wins.append(r); print('  *** WIN D=%s deficit=%d rank=%d score=%d'%(list(D),r['deficit'],r['rank'],r['score']),flush=True)
print('--- ALL %d singletons ---'%len(POOL),flush=True)
for v in POOL: rec((v,),price(make([v]),'%d'%v))
print('  done %.0fs  best gap %d'%(time.time()-t0,min(key(r)[0] for r in res)),flush=True)
print('--- ALL %d pairs ---'%(len(POOL)*(len(POOL)-1)//2),flush=True)
pairres=[]
for i,a in enumerate(POOL):
    for b in POOL[i+1:]:
        r=price(make([a,b]),'%d,%d'%(a,b)); rec((a,b),r); pairres.append((key(r),(a,b),r))
    if i%10==0: print('   %d/%d %.0fs'%(i,len(POOL),time.time()-t0),flush=True)
pairres.sort(key=lambda z:z[0])
print('  pairs done %.0fs'%(time.time()-t0),flush=True)
print('--- PRUNED triples: extensions of the best 40 pairs by (deficit-rank), plus all 3-subsets of the witness set ---',flush=True)
seen=set()
tri=[]
for _,(a,b),_ in pairres[:40]:
    for c in POOL:
        if c in (a,b): continue
        k=tuple(sorted((a,b,c)))
        if k in seen: continue
        seen.add(k); tri.append(k)
for k in itertools.combinations(WIT,3):
    k=tuple(sorted(k))
    if k not in seen: seen.add(k); tri.append(k)
print('  triples to test: %d of C(65,3)=%d  (%.1f%% coverage)'%(len(tri),65*64*63//6,100.0*len(tri)/(65*64*63//6)),flush=True)
for i,k in enumerate(tri):
    rec(k,price(make(list(k)),','.join(map(str,k))))
    if i%250==0: print('   %d/%d %.0fs'%(i,len(tri),time.time()-t0),flush=True)
json.dump(res,open('dsweep.json','w'))
res.sort(key=key)
print('\n=== BEST ROWS (sorted by deficit-rank, then score) ===',flush=True)
print('%-28s %4s %5s %5s %8s %6s %5s %8s %7s %5s'%('detach set','nz','|R|','|S|','deficit','knobs','rank','failing','score','gap'))
print('%-28s %4d %5d %5d %8d %6d %5d %8d %7d %5d  <-- WITNESS (calibration)'%(
    str(WIT),cal['nz'],cal['R'],cal['S'],cal['deficit'],cal['knobs'],cal['rank'],cal['failing'],cal['score'],CALGAP))
for r in res[:18]:
    print('%-28s %4d %5d %5d %8d %6d %5d %8d %7d %5d'%(str(r['D']),r['nz'],r['R'],r['S'],r['deficit'],r['knobs'],r['rank'],r['failing'],r['score'],r['deficit']-r['rank']))
print('\ntotal detach sets priced: %d   wins (deficit<4 or rank>7 or score>39026): %d   %.0fs'%(len(res),len(wins),time.time()-t0))
print('best score found: %d'%max(r['score'] for r in res))
