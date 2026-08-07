"""k = (structured part) + (small remainder < 2^22):  one baby table, many families by lookup."""
import ecfast as E, ec, json, time, numpy as np
p=E.p; N=E.N
P,side=ec.load(); T=ec.target()
expo=json.load(open('expo.json')); byexp={v:int(k) for k,v in expo.items()}
PT=[P[byexp[i]] for i in range(256)]
G=PT[0]
def say(s): print(s,flush=True)
def ladder(J,A,count,chunk=1<<16):
    buf=[]
    for _ in range(count):
        buf.append(J); J=E.jadd_affine(J,A)
        if len(buf)==chunk:
            for q in E.batch_norm(buf): yield q
            buf=[]
    if buf:
        for q in E.batch_norm(buf): yield q
m=1<<22
t0=time.time()
keys=np.empty(m,dtype=np.int64); i=0
for q in ladder((G[0],G[1],1),G,m):
    keys[i]=E.key(q[0]); i+=1
order=np.argsort(keys,kind='stable'); skeys=keys[order]
say('baby table 2^22 built %.0fs'%(time.time()-t0))
def lookup(Q):
    """return j in [1,2^22] with j*G == Q, else None"""
    if Q is None: return 0
    kk=E.key(Q[0]); idx=np.searchsorted(skeys,kk); out=[]
    while idx<len(skeys) and skeys[idx]==kk:
        j=int(order[idx])+1
        if E.amul(j,G)==Q: return j
        idx+=1
    return None
hits=[]
# prefix sums for runs
Pref=[None]; S=None
for i in range(256):
    S=E.aadd(S,PT[i]); Pref.append(S)
def run(a,b): return E.aadd(Pref[b],E.aneg(Pref[a]))
# F1: k = 2^i + s
t0=time.time()
for i in range(256):
    j=lookup(E.aadd(T,E.aneg(PT[i])))
    if j: hits.append(('2^i+s',i,j))
say('F1 (2^i + s<2^22) done %.0fs hits=%s'%(time.time()-t0,hits))
# F2: k = 2^i + 2^j + s
t0=time.time(); c=0
for i in range(256):
    Ti=E.aadd(T,E.aneg(PT[i]))
    for j in range(i+1,256):
        r=lookup(E.aadd(Ti,E.aneg(PT[j]))); c+=1
        if r: hits.append(('2^i+2^j+s',i,j,r))
say('F2 (2^i+2^j+s) %d lookups %.0fs hits=%s'%(c,time.time()-t0,hits))
# F3: k = run(a,b) + s
t0=time.time(); c=0
for a in range(256):
    for b in range(a+1,257):
        r=lookup(E.aadd(T,E.aneg(run(a,b)))); c+=1
        if r: hits.append(('run+s',a,b,r))
say('F3 (run + s) %d lookups %.0fs hits=%s'%(c,time.time()-t0,hits))
say('HYBRID TOTAL HITS: %s'%hits)
