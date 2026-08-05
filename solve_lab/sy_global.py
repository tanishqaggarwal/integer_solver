import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
def inv(x): return pow(x%p,-1,p)
r8=109044024666698959972204451600908701898659086097062528124234304603594878834481
r9=33371159155735472537534252650716501592825364489306217536352743247010353604716
vA=H.loadd('best_agentA_39022.json')
def base_setup():
    for v in H.freeinp: H.val[v]=vA.get(v,0)
    H.val[4287]=1
    H.val[31861]=119562606790549640390870952418684367882170154220603339634805704742270834564330392414192110
    H.val[14865]=113141528427610260107049117992526537105383080782811760722361109500341947028737388716982706
    H.val[8731]=r8; H.val[9118]=r9; H.val[9413]=0; H.val[17325]=0
    H.forward()
    H.val[4432]=H.val[19964]; H.val[7068]=H.val[2099]
    H.forward()
    H.val[950]=H.val[9106]//(13523997*p); H.val[6947]=(6122989*H.val[2239])//p; H.val[33168]=-(H.val[31731]//p)
    H.forward()
def enforce():
    # maintain G1,G2 by slaving x_4432,x_7068 after each forward
    H.forward()
    H.val[4432]=H.val[19964]+H.val[28730]
    H.val[7068]=7376877*H.val[642]+H.val[2099]
    H.forward()
base_setup(); enforce()
F=set(H.fails())
print('start fails:',len(F))
# Build closure of free inputs (2 hops) touching fails
def free_in(eqs):
    s=set()
    for i in eqs: s|=(H.eqvars[i]&H.freeinp)
    return s
clo_eqs=set(F)
for hop in range(2):
    fr=free_in(clo_eqs)
    for kn in fr:
        for i,vs in enumerate(H.eqvars):
            if kn in vs: clo_eqs.add(i)
knobs=sorted(free_in(clo_eqs) - {4432,7068})  # 4432,7068 are slaved
# exclude the pinned/used knobs we don't want to move
for x in [8731,9118]:  # keep these fixed (load residues); heal with others
    if x in knobs: knobs.remove(x)
print('closure eqs:',len(clo_eqs),' knobs:',len(knobs))
Feqs=sorted(clo_eqs)
def resids():
    enforce()
    ns={'v':H.val,'__builtins__':{}}
    return {i:eval(H.eqcode[i],ns)%p for i in Feqs}
base=resids()
nzrows=[i for i in Feqs if base[i]!=0]
print('nonzero rows in closure:',len(nzrows))
# Jacobian columns (mod p)
cols={}
for kn in knobs:
    old=H.val[kn]; H.val[kn]=old+1
    Rn=resids()
    H.val[kn]=old; enforce()
    col={i:(Rn[i]-base[i])%p for i in nzrows if (Rn[i]-base[i])%p!=0}
    if col: cols[kn]=col
knobs=[k for k in knobs if k in cols]
print('effective knobs:',len(knobs))
# Single Gaussian elimination over GF(p)
rows=nzrows[:]
Aw={i:{} for i in rows}
for kn in knobs:
    for i,v in cols[kn].items(): Aw[i][kn]=v%p
bw={i:(-base[i])%p for i in rows}
used=set(); pivot_for={}
for kn in knobs:
    prow=None
    for i in rows:
        if i in used: continue
        if Aw[i].get(kn,0)%p!=0: prow=i; break
    if prow is None: continue
    used.add(prow); pivot_for[kn]=prow
    ivp=inv(Aw[prow][kn])
    for i in rows:
        if i==prow: continue
        f=Aw[i].get(kn,0)%p
        if f==0: continue
        for kn2 in knobs:
            if kn2 in Aw[prow]:
                Aw[i][kn2]=(Aw[i].get(kn2,0)-f*ivp*Aw[prow][kn2])%p
        bw[i]=(bw[i]-f*ivp*bw[prow])%p
incon=[i for i in rows if i not in used and bw[i]%p!=0 and all(Aw[i].get(kn,0)%p==0 for kn in knobs)]
print('pivots:',len(pivot_for),' INCONSISTENT rows:',len(incon))
print('inconsistent sample:',incon[:12])
