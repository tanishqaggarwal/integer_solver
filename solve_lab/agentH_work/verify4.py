"""Close the node-triple handles over Z at the EC sum: full certification of the decompilation."""
import ec, ev, fast, json, time
from fast import St
from close3 import close
p=ec.p; S=ec.S
BITS=json.load(open('bits.json')); UW=json.load(open('uw.json'))
ALL=set(BITS['A']+BITS['B']); P,side=ec.load()
st0=St({}); U=UW['U']; W=UW['W']
FR=set(ALL)|{5096,21589,12186,16742,14853,24908,22162,30213}
res=[]
for b1,b2 in [(U[0],W[0]),(U[5],W[7]),(U[13],W[21]),(U[30],W[44])]:
    Q=ec.add(P[b1],P[b2])
    st=st0.clone().set_free({b1:1,b2:1,5096:(Q[0]-S)%p,21589:Q[1]%p})
    t0=time.time()
    out,ok,_=close(st,frozen=FR,maxsteps=400)
    nz=sorted(out.nz())
    res.append((b1,b2,out.score(),nz))
    print('b1=%d b2=%d  score=%d  nz=%s  %.1fs'%(b1,b2,out.score(),nz,time.time()-t0),flush=True)
print()
print('triple {26733,28438,32342} closed in all runs:',all(not(set(r[3])&{26733,28438,32342}) for r in res))
