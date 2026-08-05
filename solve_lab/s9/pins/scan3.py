"""Are all 512 pins satisfiable simultaneously?  And is the load map affine in the 256 bits (mod p)?"""
import sys, time, pickle, collections, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *

PINATOMS=[p_['atom'] for p_ in pins]
# key structural variables (core 1 + core 2 controls, per S9_STRUCTURE sec 3 & 9)
KEY=[12186,14853,22162,16742,24908,30213,24453,7068,4432,2099,19964,
     29322,3558,35389,6671,11150,25739,37758,15298,
     30454,10261,16787,25199,21344,21589,5096,
     18123,17576,25614,34220,21202,32453,15286,38170,8599,21839]

def act(v, bits, pinsclosed=True):
    s={}
    for b in BITS:
        for pn in bitpins[b]: s[pn['h']]=0
    for b in bits:
        s[b]=1
        if pinsclosed:
            for pn in bitpins[b]: s[pn['B']]=pn['HUGE']
    ripple(v,s)

def snap(v): return {k: v[k] % P for k in KEY}

if __name__=='__main__':
    t0=time.time()
    v0=list(BASE); act(v0, [2081,24601])
    bad=[a for a in PINATOMS if evalpoly(polys[a],v0)!=0]
    print('pins violated at base state:', len(bad))
    base=snap(v0)

    # all-bits-on: are all 512 pins simultaneously satisfiable?
    vall=list(BASE); act(vall, BITS)
    badall=[a for a in PINATOMS if evalpoly(polys[a],vall)!=0]
    print('pins violated with ALL 256 bits on + all pins closed:', len(badall))
    print('   (=> pin set is jointly satisfiable as a constraint family)' if not badall else badall[:10])

    # per-bit deltas of the key variables mod p
    deltas={}
    for b in BITS:
        if BASE[b]: continue
        v=list(BASE); act(v,[2081,24601,b])
        s=snap(v)
        deltas[b]={k:(s[k]-base[k])%P for k in KEY if (s[k]-base[k])%P}
    pickle.dump((base,deltas), open('pins/deltas.pkl','wb'))
    print(f'\nper-bit key-variable deltas computed [{time.time()-t0:.0f}s]')
    cnt=collections.Counter()
    for b,d in deltas.items():
        for k in d: cnt[k]+=1
    print('how many bits move each key var (mod p):')
    for k in KEY: print(f'   x_{k}: {cnt[k]}')

    # additivity test on random pairs
    random.seed(0); cand=[b for b in BITS if not BASE[b]]
    print('\nadditivity test (mod p) on random bit pairs:')
    for _ in range(8):
        b1,b2=random.sample(cand,2)
        v=list(BASE); act(v,[2081,24601,b1,b2]); s=snap(v)
        ok=True; diffs=[]
        for k in KEY:
            pred=(base[k]+deltas[b1].get(k,0)+deltas[b2].get(k,0))%P
            if pred!=s[k]: ok=False; diffs.append(k)
        print(f'   ({b1},{b2}): affine={ok} {diffs[:6]}')
