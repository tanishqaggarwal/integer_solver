"""Structured bit-pattern families for k (each is a subset of exponents 0..255)."""
import ecfast as E, ec, json, time
p=E.p; N=E.N
P,side=ec.load(); T=ec.target()
expo=json.load(open('expo.json')); byexp={v:int(k) for k,v in expo.items()}
PT=[P[byexp[i]] for i in range(256)]
# prefix sums of powers: Pref[i] = P_0+..+P_{i-1}
Pref=[None]
S=None
for i in range(256):
    S=E.aadd(S,PT[i]); Pref.append(S)
def sumset(idxs):
    R=None
    for i in idxs: R=E.aadd(R,PT[i])
    return R
hits=[]
t0=time.time()
# 1. periodic patterns: bit i set iff (i mod d) in R, for d<=12 and all nonempty R subsets, window [lo,hi)
tested=0
for d in range(1,13):
    for mask in range(1,1<<d):
        Rs=[r for r in range(d) if mask>>r&1]
        for lo in range(0,256,8):
            for hi in (128,192,256):
                if hi<=lo: continue
                idxs=[i for i in range(lo,hi) if (i%d) in Rs]
                if not idxs or len(idxs)>256: continue
                q=sumset(idxs); tested+=1
                if q==T: hits.append(('periodic',d,mask,lo,hi))
    if d<=8: print('  d=%d tested=%d %.0fs'%(d,tested,time.time()-t0),flush=True)
print('periodic families tested:',tested,'hits',hits,'%.1fs'%(time.time()-t0))
# 2. two runs: [a,b) U [c,e)
t0=time.time(); cnt=0
for a in range(0,256,2):
    for b in range(a+1,257,2):
        base=E.aadd(Pref[b],E.aneg(Pref[a]))
        for c in range(b,256,4):
            for e in range(c+1,257,4):
                q=E.aadd(base,E.aadd(Pref[e],E.aneg(Pref[c]))); cnt+=1
                if q==T: hits.append(('tworuns',a,b,c,e))
    if a%64==0: print('  tworuns a=%d cnt=%d %.0fs'%(a,cnt,time.time()-t0),flush=True)
print('two-run families tested:',cnt,'hits',hits,'%.1fs'%(time.time()-t0))
print('TOTAL HITS:',hits)
