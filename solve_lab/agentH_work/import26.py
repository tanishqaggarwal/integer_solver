"""Load the 39,026 witness, express it in MY frame, and hill-climb with my extra free inputs."""
import ev, fast, json, os, time, sys
from fast import St, chk
from collections import defaultdict
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
W=json.load(open('../best/new_instance_partial_39026.json'))
v=[0]*38748
for k,val in W.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
FREE=ev.F['free0']
fv={u:v[u] for u in FREE if v[u]!=0}
print('witness free-input values in my frame: %d nonzero of %d'%(len(fv),len(FREE)))
st=St(fv)
print('forward-eval of the witness free inputs in MY frame: score',st.score(),'nz atoms',len(st.nz()))
print('  nz:',sorted(st.nz())[:20])
print('  failing:',sorted(st.fails)[:20])
# how far is my reconstruction from the witness itself?
diff=[i for i in range(38748) if st.v[i]!=v[i]]
print('variables where my forward eval differs from the witness:',len(diff),diff[:20])
json.dump({str(a):b for a,b in fv.items()},open('w26_fv.json','w'))
