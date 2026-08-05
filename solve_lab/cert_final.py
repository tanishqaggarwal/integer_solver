import heal_harness as H
import json, pickle
p=H.p
SCR='/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad'
ATOMS=[]; reprs=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        dd=json.loads(line); ATOMS.append([(tuple(m),c) for m,c in dd['poly']]); reprs.append(dd.get('repr',''))
cq=pickle.load(open(SCR+'/cert_quadrant.pkl','rb')); cert=cq['cert']
d=H.loadd('best/new_instance_partial_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
print(f"Tightest quadrant certificate: {len(cert)} atoms (control bits fixed, boolean atoms excluded)")
print("Atoms and their role:")
# classify: gap, load, pin (control bit * (target-CONST)), other
pins=[]
for a,mv in sorted(cert):
    rp=reprs[a]
    role=""
    if a in (20862,20864): role="GAP (G1/G2 target)"
    elif a in (18081,18084,29377,35321): role="LOAD atom"
    elif '*' in rp and ' - ' in rp and any(str(x) in rp for x in range(10)):
        # detect pin: var*(var - HUGECONST)
        role="pin/product"
    print(f"  atom {a:5d} [{role:18s}]: {rp[:78]}")
# identify control bits appearing as selectors in cert pin atoms
print("\nControl-bit selectors in certificate pin atoms (flip candidates for quadrant search):")
ctrlbits=set()
for a,mv in cert:
    for m,c in ATOMS[a]:
        if len(m)==2:
            # bilinear a*b: check if one is a boolean control bit (value 0/1)
            for v in m:
                if H.val[v] in (0,1) and v!=15298:
                    ctrlbits.add(v)
# Also known override bits
for b in [2081,4287,24601,13195]:
    if any(b in set(x for mm,cc in ATOMS[a] for x in mm) for a,_ in cert): ctrlbits.add(b)
print("  candidate control bits:", sorted(ctrlbits))
for b in sorted(ctrlbits):
    print(f"    x_{b} = {H.val[b]}  free={b in H.freeinp}")

# EXACT verification of the certificate (recompute Jacobian rows and check LHS=0, RHS!=0)
print("\n=== exact re-verification of certificate ===")
from collections import defaultdict
gatepoly=pickle.load(open(SCR+'/gatepoly.pkl','rb')); gate_items={t:list(pl.items()) for t,pl in gatepoly.items()}
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
def boolean_var(i):
    vs=set(x for m,c in ATOMS[i] for x in m)
    if len(vs)!=1: return None
    v=next(iter(vs)); c1=c2=0
    for m,c in ATOMS[i]:
        if len(m)==2: c2+=c
        elif len(m)==1: c1+=c
        elif len(m)==0 and c!=0: return None
    if c2 and c1 and (c1+c2)%p==0: return v
    return None
BOOLFREE=set(bv for i in range(len(ATOMS)) if (bv:=boolean_var(i)) is not None)&H.freeinp
val=H.val
def ppi(items,v):
    s=0
    for m,c in items:
        cnt=m.count(v)
        if cnt==0: continue
        term=(c*cnt)%p; seen=False
        for u in m:
            if u==v and not seen: seen=True; continue
            term=term*val[u]%p
        s=(s+term)%p
    return s
def av(i):
    s=0
    for m,c in ATOMS[i]:
        tt=c%p
        for v in m: tt=tt*val[v]%p
        s=(s+tt)%p
    return s
localJac={}
for t,items in gate_items.items():
    localJac[t]={v:ppi(items,v) for v in set(x for m,c in items for x in m)}
KNOBS=set(f for f in H.freeinp if f not in BOOLFREE)
def jrow(a):
    # Jacobian row of atom a wrt KNOBS via reverse? use forward per needed... do full forward per knob is heavy;
    # instead compute directly: dA/dknob = sum over vars dA/dvar * dvar/dknob. Use one forward-mode per atom's cone
    # Simpler: reuse global forward-mode but only accumulate this atom.
    jac={}
    for v in set(x for m,c in ATOMS[a] for x in m):
        pv=ppi(ATOMS[a],v)
        if pv: jac[v]=pv
    row=defaultdict(int)
    # forward-mode from each knob is expensive; do reverse via adjoint over gates in cone
    # adjoint approach:
    adj=dict(jac)
    for k in range(len(H.order)-1,-1,-1):
        t=H.order[k]
        at=adj.get(t)
        if not at: continue
        for v,cf in localJac[t].items():
            adj[v]=(adj.get(v,0)+at*cf)%p
    for f in KNOBS:
        if adj.get(f): row[f]=adj[f]%p
    return dict(row)
lhs=defaultdict(int); rhs=0
for a,mv in cert:
    row=jrow(a)
    ra=0 if a not in (20862,20864) else av(a)  # gap rows carry residual on RHS side
    for c,cf in row.items(): lhs[c]=(lhs[c]+mv*cf)%p
    # RHS: for gap target rows the system rhs is -residual; for kept rows it's 0
    if a in (20862,20864): rhs=(rhs+mv*((-av(a))%p))%p
lhs={c:v for c,v in lhs.items() if v%p}
print("certificate LHS nonzero knob-columns:", len(lhs), " (should be 0)")
print("certificate RHS (should be != 0):", rhs%p!=0, "value:", rhs%p)
