import pickle, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
NV=38748
boolv=set(pickle.load(open('boolvars.pkl','rb')))
freeinp=[x for x in range(NV) if x not in definer]
bfree=[b for b in freeinp if b in boolv]
v0=H.load_assignment('../best/new_instance_partial_39022.json')
hits=[]
for b in bfree:
    v=list(v0); ripple(v,{b:1-v0[b]})
    if v[8599]==1 and v[21839]==1:
        hits.append((b, v[38170], v[15298], v[2754], v[24673]))
print(f'free bits that give x_8599=1 while keeping x_21839=1: {len(hits)}')
for h in hits[:25]: print('   x_%d -> x_38170=%d x_15298=%d x_2754=%d x_24673=%d'%h)
pickle.dump([h[0] for h in hits], open('hits8599.pkl','wb'))
