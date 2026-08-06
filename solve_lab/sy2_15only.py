import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
def inv(x): return pow(x%p,-1,p)
fc=H.loadd('fc_partial.json')
for v in H.freeinp: H.val[v]=fc.get(v,0)
H.val[14865]=H.val[12553]; H.val[31861]=H.val[6418]; H.val[33168]=0; H.val[950]=0; H.val[6947]=0
H.forward()
if H.val[37720]%9994531==0: H.val[8976]=H.val[37720]//9994531
H.forward()
cons=[287,1531,3081,7425,8273,8470,9708,13790,18030,27706,29926,30383,32138,37297,38740]
# fine-grained knobs = free inputs in cons whose perturbation changes SOME eq mod p
knobs=set()
for i in cons: knobs|=(H.eqvars[i]&H.freeinp)
knobs-={4287,2081,24601,6418,12553,8731,9118,4432,7068}  # protect loads/G1/G2/x_6418/x_12553; allow x_14865,x_31861
knobs=sorted(knobs)
def rc():
    ns={'v':H.val,'__builtins__':{}}
    return {i:eval(H.eqcode[i],ns)%p for i in cons}
base=rc()
nz=[i for i in cons if base[i]!=0]
print('nonzero consumers mod p:',len(nz))
# fine-grained knobs only (nonzero mod-p column on cons)
cols={}; fine=[]
for kn in knobs:
    H.val[kn]+=1; H.forward(); Rn=rc(); H.val[kn]-=1; H.forward()
    col={i:(Rn[i]-base[i])%p for i in cons if (Rn[i]-base[i])%p!=0}
    if col: cols[kn]=col; fine.append(kn)
print('fine-grained knobs affecting consumers mod p:',len(fine),fine[:25])
# Gaussian over the 15 cons rows only
rowdata={i:{} for i in cons}
for kn,col in cols.items():
    for i,v in col.items(): rowdata[i][kn]=v
rhs={i:(-base[i])%p for i in cons}
used=set(); piv=[]
for kn in fine:
    prow=None
    for i in cons:
        if i in used: continue
        if rowdata[i].get(kn,0)%p!=0: prow=i;break
    if prow is None: continue
    used.add(prow); piv.append((kn,prow))
    ipv=inv(rowdata[prow][kn])
    for c in list(rowdata[prow]): rowdata[prow][c]=(rowdata[prow][c]*ipv)%p
    rhs[prow]=(rhs[prow]*ipv)%p
    for i in cons:
        if i==prow: continue
        f=rowdata[i].get(kn,0)%p
        if f==0: continue
        for c,val in rowdata[prow].items(): rowdata[i][c]=(rowdata[i].get(c,0)-f*val)%p
        rhs[i]=(rhs[i]-f*rhs[prow])%p
incon=[i for i in cons if i not in used and rhs[i]%p!=0]
print('rank over 15 consumers:',len(piv),' inconsistent (of 15):',len(incon), incon)
