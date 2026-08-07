"""Decisive test: does setting X3,Y3 := (P_b1 + P_b2) make the chord residuals A,B vanish mod p?"""
import ec, ev, fast, json, time
from fast import St
from close2 import close
p=ec.p; S=ec.S
BITS=json.load(open('bits.json')); UW=json.load(open('uw.json'))
ALL=set(BITS['A']+BITS['B']); P,side=ec.load()
st0=St({}); U=UW['U']; W=UW['W']
FR=set(ALL)|{5096,21589,18956,24468,7497,11436,22820,14393,12186,16742,14853,24908,22162,30213}
tests=[(U[i],W[j]) for i in range(0,len(U),11) for j in range(0,len(W),13)]
good=0; ctrl=0
for b1,b2 in tests:
    Q=ec.add(P[b1],P[b2])
    for tag,dx in (('exact',0),('off-by-1',1)):
        st=st0.clone().set_free({b1:1,b2:1,5096:(Q[0]-S+dx)%p,21589:Q[1]%p})
        out,_,_=close(st,frozen=FR,maxsteps=400)
        A=out.v[25614]%p; B=out.v[34220]%p        # chord residuals
        z=(A==0 and B==0)
        if tag=='exact' and z: good+=1
        if tag=='off-by-1' and not z: ctrl+=1
print('pairs tested:',len(tests))
print('  A==B==0 mod p with X3,Y3 := P_b1+P_b2 :', good,'/',len(tests))
print('  A,B nonzero when X3 perturbed by 1     :', ctrl,'/',len(tests))
