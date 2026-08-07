"""Priority 2: certify the EC decompilation.  Two bits in DIFFERENT root groups, X3/Y3 UNFROZEN,
set to the predicted EC sum; check every check-atom of the node triple closes."""
import ec, ev, fast, json, sys, time
from fast import St
from close2 import close
p=ec.p; S=ec.S
BITS=json.load(open('bits.json')); UW=json.load(open('uw.json'))
ALL=set(BITS['A']+BITS['B'])
P,side=ec.load()
st0=St({})
U=UW['U']; W=UW['W']
ok=bad=0
rows=[]
t0=time.time()
tests=[(U[i],W[j]) for i in range(0,len(U),7) for j in range(0,len(W),9)]
for b1,b2 in tests:
    Q=ec.add(P[b1],P[b2])            # predicted delivered point
    X3=(Q[0]-S)%p; Y3=Q[1]%p         # un-shift back to instance coordinates
    st=st0.clone().set_free({b1:1,b2:1,5096:X3,21589:Y3})
    out,_,_=close(st,frozen=set(ALL)|{5096,21589,18956,24468,7497,11436,22820,14393,
                                      12186,16742,14853,24908,22162,30213},maxsteps=400)
    nz=set(out.nz())
    triple_closed = not (nz & {26733,28438,32342})
    got=((out.v[23927]+S)%p, out.v[19083]%p)
    rows.append((b1,b2,triple_closed,got==Q,sorted(nz)))
    if triple_closed and got==Q: ok+=1
    else: bad+=1
print('tests %d  certified %d  failed %d   %.1fs'%(len(tests),ok,bad,time.time()-t0))
for r in rows[:6]: print('  b1=%d b2=%d tripleclosed=%s delivered==P1+P2:%s nz=%s'%r)
if bad:
    for r in rows:
        if not (r[2] and r[3]): print('  FAIL',r); break
# control: a WRONG X3/Y3 must leave the triple open
b1,b2=tests[0]
Q=ec.add(P[b1],P[b2])
st=st0.clone().set_free({b1:1,b2:1,5096:(Q[0]-S+1)%p,21589:Q[1]})
out,_,_=close(st,frozen=set(ALL)|{5096,21589,18956,24468,7497,11436,22820,14393,
                                  12186,16742,14853,24908,22162,30213},maxsteps=400)
print('CONTROL (X3 off by 1): triple atoms nonzero =', sorted(set(out.nz())&{26733,28438,32342}))
