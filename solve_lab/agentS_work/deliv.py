import sys, json, collections
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E
P=C.P
FOOT=collections.defaultdict(set)
for e,(issq,outer,terms) in enumerate(H.eqt):
    for c,a in terms:
        if a>=0: FOOT[a].add(e)
NF={a:len(s) for a,s in FOOT.items()}
d=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v=[0]*E.NV
for k,val in d.items(): v[int(k.split('_')[1])]=int(val)
bad=E.badatoms(v); ff=E.eqfails(bad)
print("deliverable: %d nonzero atoms, %d failing eqs -> score %d"%(len(bad),len(ff),39033-len(ff)))
print("failing eqs:",ff)
for a in sorted(bad):
    print(f"  a{a}: nf={NF[a]} resid_bits={bad[a].bit_length()} mod p={'0' if bad[a]%P==0 else 'nonzero'}  eqs={sorted(FOOT[a])}")
print("\nunion of footprints:",sorted(set().union(*[FOOT[a] for a in bad])))
# how many footprint-1 atoms exist
f1=[a for a in NF if NF[a]==1]
print("\natoms with footprint 1:",len(f1))
# do footprint-1 atoms cover distinct equations?
eqs1=collections.Counter()
for a in f1: eqs1[next(iter(FOOT[a]))]+=1
print("distinct equations covered by footprint-1 atoms:",len(eqs1))
