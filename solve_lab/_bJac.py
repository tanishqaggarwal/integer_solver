import heal_harness as H, re, time
from collections import defaultdict
p=H.p
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
base[4287]=1; base[2081]=1; base[9118]=base[7068]; base[8731]=base[4432]
for v in H.freeinp: H.val[v]=base[v]
H.forward()
F=H.fails()
print(f"branch B fails: {len(F)}")
# free vars appearing in failing eqs
Kset=sorted(set(v for i in F for v in H.eqvars[i] if v in H.freeinp))
print(f"free vars in failing eqs: {len(Kset)}")
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
touch={v:sorted(set(desc_of[v])) for v in Kset}
def setf(v,x):
    H.val[v]=x
    for k in touch[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
def residp():
    return [eval(H.eqcode[i],ns)%p for i in F]
b0=residp()
inv2=pow(2,p-2,p)
J=[[0]*len(Kset) for _ in F]
for j,v in enumerate(Kset):
    o=base[v]; setf(v,o+1); rp=residp(); setf(v,o-1); rm=residp(); setf(v,o)
    for k in range(len(F)): J[k][j]=((rp[k]-rm[k])*inv2)%p
def rank_modp(M):
    M=[row[:] for row in M]; rows=len(M); cols=len(M[0]) if rows else 0; r=0
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if M[i][c]%p: piv=i;break
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; iv=pow(M[r][c]%p,p-2,p)
        M[r]=[(x*iv)%p for x in M[r]]
        for i in range(rows):
            if i!=r and M[i][c]%p:
                f=M[i][c]; M[i]=[(M[i][t]-f*M[r][t])%p for t in range(cols)]
        r+=1
    return r
rJ=rank_modp(J); rJb=rank_modp([J[k]+[(-b0[k])%p] for k in range(len(F))])
print(f"mod-p Jacobian: rank(J)={rJ}, rank([J|-b])={rJb} -> {'CONSISTENT (Newton mod p exists!)' if rJ==rJb else 'INCONSISTENT'}")
print(f"num free knobs={len(Kset)}, num eqs={len(F)}")
