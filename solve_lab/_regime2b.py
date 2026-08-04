import heal_harness as H, json, pickle
p=H.p
C=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/atomcache.pkl','rb'))
atoms=C['atoms']; eq2atoms=C['eq2atoms']
def av(ai):
    a=atoms[ai]; s=0
    for vl,c in a['poly']:
        t=c
        for v in vl: t*=H.val[v]
        s+=t
    return s
d2=H.loadd('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/regime2.json')
def setfree(dd):
    for v in range(H.NVARS): H.val[v]=dd.get(v,0)
setfree(d2); H.forward()
# handle vars: are they free?
for v in [30317,2936,5146,5101,32017,26789]:
    print(f"x_{v}: free={v in H.freeinp} val={H.val[v]}")
# L values (mod-p reduced = 0, get full integers)
L1=H.val[11150]; L2=H.val[25739]; L3=H.val[37758]
print(f"\nL1 % p == 0: {L1%p==0};  L1/p int? {L1%p==0}")
print(f"L2 % p == 0: {L2%p==0}")
print(f"L3 % p == 0: {L3%p==0}")
print(f"L2/p mod 6672769 = {(L2//p)%6672769}")
# Set handles
d3=dict(d2)
d3[30317]=(-L1//p)      # x_4007 = p*x_30317 = -L1
d3[2936]=(537773*L3//p) # x_35605 = p*x_2936 = 537773*L3
# x_5146 only if L2 % (6672769*p)==0
m=6672769
if L2 % (m*p)==0:
    d3[5146]=L2//(m*p)
    print("L2 divisible by 6672769*p, x_5146 set")
else:
    print(f"L2 NOT divisible by 6672769*p (need L2/p % {m} = 0, got {(L2//p)%m})")
    d3[5146]=d2.get(5146,0)
setfree(d3); H.forward()
print(f"\n  x_4007={H.val[4007]}  should be -L1={-L1}  match={H.val[4007]==-L1}")
print(f"  x_35605={H.val[35605]}  should be 537773*L3={537773*L3}  match={H.val[35605]==537773*L3}")
F=sorted(H.fails())
print(f"  fails: {len(F)}: {F}")
# census
from collections import defaultdict
atom_eqs=defaultdict(list)
for e in F:
    for ai in eq2atoms.get(e,[]):
        if av(ai)!=0: atom_eqs[ai].append(e)
print(f"\n  nonzero atoms: {len(atom_eqs)}")
for ai in sorted(atom_eqs):
    print(f"    atom#{ai} n_eq={atoms[ai]['n_eq']} fails={atom_eqs[ai][:8]}: {atoms[ai]['repr'][:75]}")
json.dump({f"x_{k}":str(v) for k,v in d3.items()},open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/regime3.json','w'))
