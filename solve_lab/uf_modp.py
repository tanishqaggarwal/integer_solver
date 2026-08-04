import json
import heal_harness as H
p=H.p
# weighted union-find mod p: x_a ≡ mult*x_root + off
parent=list(range(H.NVARS)); mult=[1]*H.NVARS; off=[0]*H.NVARS
def find(x):
    if parent[x]==x: return x,1,0
    r,m,o=find(parent[x])
    parent[x]=r; mult[x]=(mult[x]*m)%p; off[x]=(off[x]*m+o)%p
    return r,mult[x],off[x]
contradictions=[]
def union(a,b,m_ab,o_ab):  # x_a ≡ m_ab*x_b + o_ab mod p
    ra,ma,oa=find(a); rb,mb,ob=find(b)
    # x_a ≡ ma*ra+oa ; x_b ≡ mb*rb+ob ; want x_a=m_ab*x_b+o_ab
    # ma*ra+oa = m_ab*(mb*rb+ob)+o_ab
    if ra==rb:
        # check consistency: ma*ra+oa must equal m_ab*(mb*ra+ob)+o_ab for the SAME root ra
        lhs_m=ma; lhs_o=oa; rhs_m=(m_ab*mb)%p; rhs_o=(m_ab*ob+o_ab)%p
        if lhs_m!=rhs_m or lhs_o!=rhs_o:
            contradictions.append((a,b))
        return
    # attach rb under ra: ra ≡ (x_a-oa)/ma ; express rb via ra
    # x_a=m_ab*x_b+o_ab -> ma*ra+oa=m_ab*mb*rb+m_ab*ob+o_ab -> rb=(ma*ra+oa-m_ab*ob-o_ab)/(m_ab*mb)
    inv=pow((m_ab*mb)%p,-1,p)
    parent[rb]=ra; mult[rb]=(ma*inv)%p; off[rb]=((oa-m_ab*ob-o_ab)*inv)%p
# load atoms + current quadrant bit values
d=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
# wire members (=p): value % p == 0 and it's the union-find class... use: var is wire if H.val says ~p. 
# Simpler: treat a product term x_w*x_p as p-multiple if H.val[x_w]%p==0
used=0
for a in atoms:
    poly=a['poly']
    # 2-term linear: [[x_a],c1],[[x_b],c2]  or with const
    terms=[(tuple(vs),c) for vs,c in poly]
    lin=[(vs[0],c) for vs,c in terms if len(vs)==1]
    const=sum(c for vs,c in terms if len(vs)==0)
    quad=[(vs,c) for vs,c in terms if len(vs)==2]
    # case A: pure 2-term linear a*x + b*y + const =0 (no quad) -> x ≡ ... 
    if not quad and len(lin)==2 and (lin[0][1]%p) and (lin[1][1]%p):
        (va,ca),(vb,cb)=lin
        # ca*va+cb*vb+const=0 -> va ≡ (-cb/ca)vb + (-const/ca)
        m=(-cb*pow(ca,-1,p))%p; o=(-const*pow(ca,-1,p))%p
        union(va,vb,m,o); used+=1
    # case B: 1-term linear + const = pin/lock: ca*va+const=0 -> va ≡ -const/ca
    elif not quad and len(lin)==1 and (lin[0][1]%p):
        va,ca=lin[0]; val=(-const*pow(ca,-1,p))%p
        # union with a virtual constant node? use var itself pinned: store as union to a fixed sentinel
        ra,ma,oa=find(va)
        # ma*ra+oa ≡ val -> ra ≡ (val-oa)/ma
        # pin root ra to constant: emulate by unioning to node 0 with offset? skip—track pins separately
        pass
    # case C: scaled diff minus p-partner: ca*va + cb*vb + (quad terms that are p-multiples)+const=0
    elif quad and len(lin)==2:
        # check all quad terms are p-multiples at current val
        allp=all(H.val[vs[0]]%p==0 or H.val[vs[1]]%p==0 for vs,c in quad)
        if allp and (lin[0][1]%p) and (lin[1][1]%p) and const%p==0:
            (va,ca),(vb,cb)=lin
            m=(-cb*pow(ca,-1,p))%p; o=0
            union(va,vb,m,o); used+=1
print(f"used {used} mod-p linear relations")
# check core vars
for a,b in [(14853,12186),(16742,24908)]:
    ra,ma,oa=find(a); rb,mb,ob=find(b)
    if ra==rb:
        # x_a ≡ ma*r+oa, x_b ≡ mb*r+ob ; are they forced equal? need ma=mb,oa=ob for all r? 
        # x_a - x_b ≡ (ma-mb)*r+(oa-ob). Forced equal iff ma==mb and oa==ob
        forced_equal = (ma==mb and oa==ob)
        print(f"x_{a},x_{b}: SAME class. forced x_a-x_b ≡ ({(ma-mb)%p})*root+{(oa-ob)%p}. Can be equal? {'NO-PINNED-DIFFERENT' if (ma==mb and oa!=ob) else 'depends on root' if ma!=mb else 'YES(already equal)'}")
    else:
        print(f"x_{a},x_{b}: DIFFERENT classes -> NOT forced, core condition FREE to satisfy")
print(f"contradictions found: {len(contradictions)}")
