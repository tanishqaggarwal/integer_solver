import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import frameB
t0=time.time()
fr = frameB.Frame([642,28730,29854,31864])
print('frame built %.1fs free=%d checks=%d' % (time.time()-t0, len(fr.free), len(fr.checks)))
W = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v=[0]*frameB.NV
for k,val in W.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
fv={u:v[u] for u in fr.free if v[u]!=0}
st=frameB.State(fr,fv)
print('score',st.score(),'failing',sorted(st.fails))
print('nz atoms',sorted(st.nz()))
dif=[i for i in range(frameB.NV) if st.v[i]!=v[i]]
print('vars differing:',len(dif))
