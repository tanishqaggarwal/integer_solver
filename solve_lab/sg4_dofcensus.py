import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
V0=H.val[:]
baseF=set(H.fails())
def R1R2():
    return (H.val[2099]-H.val[7068])%p, (H.val[19964]-H.val[4432])%p
r1_0,r2_0=R1R2()
print(f"baseline R1={r1_0}  R2={r2_0}  ({len(baseF)} fails)")
# candidate free inputs affecting R1,R2 + broaden to all free ancestors of x_2099,x_19964,x_7068,x_4432
cands=set([2081,4287,6418,9118,31861,8731,12553,14865,7068,4432])
for t in [2099,19964,7068,4432,642,28730,17601,24908]:
    cands|= {a for a in H.anc.get(t,set()) if a in H.freeinp}
    if t in H.freeinp: cands.add(t)
print(f"candidate free inputs: {len(cands)}")
print(f"{'var':>8} {'dR1modp':>10} {'dR2modp':>10} {'fanout':>7}  note")
rows=[]
for v in sorted(cands):
    # perturb by 1
    H.val[v]=V0[v]+1
    H.forward()
    dr1=(H.val[2099]-H.val[7068]-r1_0-(0))
    r1n,r2n=R1R2()
    dR1=(r1n-r1_0)%p; dR2=(r2n-r2_0)%p
    newF=set(H.fails())
    fo=len(newF-baseF)
    rows.append((fo,v,dR1,dR2))
    # restore
    for k in range(len(H.val)): H.val[k]=V0[k]
    H.forward()
rows.sort()
for fo,v,dR1,dR2 in rows:
    s1 = '0' if dR1==0 else ('nz' if dR1!=0 else '?')
    note=''
    if dR1!=0 or dR2!=0:
        note='SHIFTS-RESIDUE'
    if fo==0: note+=' [0-FANOUT]'
    print(f"x_{v:>6} {('0' if dR1==0 else 'NONZERO'):>10} {('0' if dR2==0 else 'NONZERO'):>10} {fo:>7}  {note}")
# Summary: any 0-fanout knob that shifts a residue?
clean=[(v,dR1,dR2) for fo,v,dR1,dR2 in rows if fo==0 and (dR1 or dR2)]
print(f"\n0-FANOUT residue-shifting knobs: {clean}")
