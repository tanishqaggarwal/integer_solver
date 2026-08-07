"""Can the 39,026 deliverable itself be improved by pure-handle repair?"""
import sys, json, collections, time
sys.path.insert(0,'.')
import common as C, hrepair as HR
import harness as H, engine as E
P=C.P
d=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vD=[0]*E.NV
for k,val in d.items(): vD[int(k.split('_')[1])]=int(val)
badD=E.badatoms(vD); ffD=E.eqfails(badD)
print("deliverable as given: %d fails, score %d, atoms %s"%(len(ffD),39033-len(ffD),sorted(badD)))
# re-derive from its free-variable values through my forward map
seed={f:vD[f] for f in sorted(E.FREE)}
t0=time.time(); v2=E.forward(seed); bad2=E.badatoms(v2); ff2=E.eqfails(bad2)
same=sum(1 for i in range(E.NV) if v2[i]==vD[i])
print("forward(free of deliverable): %d/%d vars identical, %d fails, score %d, atoms %s (%.0fs)"
      %(same,E.NV,len(ff2),39033-len(ff2),sorted(bad2),time.time()-t0))
# handle repair on that
n,bad3,v3,ns=HR.repair(seed)
print("after pure-handle repair: %d fails, score %d, atoms %s"%(n,39033-n,sorted(bad3)))
for a in sorted(bad3):
    h=HR.find_handle(a)
    print(f"   a{a} nf={HR.NF[a]} handle={h[0] if h else None} step={'p' if h and abs(h[1])==P else (h[1] if h else None)} R%p={'0' if bad3[a]%P==0 else 'nz'}")
HR.save()
if n<7:
    json.dump({f"x_{i}":int(v3[i]) for i in range(E.NV) if v3[i]!=0}, open('S_improved_%d.json'%(39033-n),'w'))
    print("WROTE S_improved_%d.json"%(39033-n))
