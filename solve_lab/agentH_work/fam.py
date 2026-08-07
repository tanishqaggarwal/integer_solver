"""Remaining structured-k families: two-run exhaustive, periodic patterns, BSGS on -T."""
import ecfast as E, ec, json, time, numpy as np, sys
p=E.p; N=E.N
P,side=ec.load(); T=ec.target()
expo=json.load(open('expo.json')); byexp={v:int(k) for k,v in expo.items()}
PT=[P[byexp[i]] for i in range(256)]
def say(s): print(s,flush=True)

# prefix sums
Pref=[None]; S=None
for i in range(256):
    S=E.aadd(S,PT[i]); Pref.append(S)
def run(a,b):  # sum of exponents [a,b)
    return E.aadd(Pref[b],E.aneg(Pref[a]))

# ---- family 1: exactly two runs of consecutive set bits (exhaustive) ----
t0=time.time()
tab={}
for a in range(257):
    for b in range(a+1,257):
        q=run(a,b)
        if q is not None: tab.setdefault(q,(a,b))
say('one-run table: %d distinct sums  %.1fs'%(len(tab),time.time()-t0))
hits=[]
if T in tab: hits.append(('onerun',tab[T]))
for q,(a,b) in tab.items():
    r=E.aadd(T,E.aneg(q))
    if r in tab:
        c,d=tab[r]
        if b<=c or d<=a: hits.append(('tworun',(a,b),(c,d)))
say('TWO-RUN family (exhaustive over all 33k x 33k): hits = %s   %.1fs'%(hits[:3],time.time()-t0))

# ---- family 2: periodic bit patterns, period d <= 16, all residue sets ----
t0=time.time(); phits=[]
for d in range(1,17):
    Q=[]
    for r in range(d):
        s=None
        for j in range(r,256,d): s=E.aadd(s,PT[j])
        Q.append(s)
    # gray-code over subsets of residues
    cur=None; prev=0
    for mask in range(1,1<<d):
        g=mask^(mask>>1)
        diff=g^prev
        b=diff.bit_length()-1
        cur=E.aadd(cur,Q[b]) if (g>>b&1) else E.aadd(cur,E.aneg(Q[b]))
        prev=g
        if cur==T: phits.append(('periodic',d,g))
    say('  d=%d done  %.1fs'%(d,time.time()-t0))
say('PERIODIC family (d<=16, all 2^d residue sets): hits = %s   %.1fs'%(phits,time.time()-t0))

# ---- family 3: k = N - c with c < 2^44  (BSGS on -T) ----
import dl as DL
say('--- BSGS on -T : covers k = N - c, c < 2^44 ---')
Tn=E.aneg(T)
m=1<<22
keys=np.empty(m,dtype=np.int64)
G=PT[0]
J=(G[0],G[1],1); i=0
def ladder(start_jac, add_affine, count, chunk=1<<16):
    Jl=start_jac; buf=[]
    for _ in range(count):
        buf.append(Jl); Jl=E.jadd_affine(Jl,add_affine)
        if len(buf)==chunk:
            for q in E.batch_norm(buf): yield q
            buf=[]
    if buf:
        for q in E.batch_norm(buf): yield q
t0=time.time()
for q in ladder(J,G,m):
    keys[i]=E.key(q[0]); i+=1
order=np.argsort(keys,kind='stable'); skeys=keys[order]
say('  baby built+sorted %.0fs'%(time.time()-t0))
mG=E.amul(m,G); nmG=E.aneg(mG)
found=None; gi=0
for q in ladder((Tn[0],Tn[1],1),nmG,1<<22):
    kk=E.key(q[0]); idx=np.searchsorted(skeys,kk)
    while idx<len(skeys) and skeys[idx]==kk:
        j=int(order[idx])+1
        cand=(gi*m+j)%N
        if E.amul(cand,G)==Tn: found=cand; break
        idx+=1
    if found is not None: break
    gi+=1
say('NEG-BSGS: c = %s  (k = N-c)  %.0fs'%(found,time.time()-t0))
if found is not None:
    k=(N-found)%N
    say('*** k = %d ***'%k)
