"""Neighborhood of the actual deliverable: keep its 30 non-leaf free settings, swap its leaf pair
   {24601 (a-side), 2081 (b-side)} for other pairs, score with E.forward."""
import sys,json,time,pickle,random,collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentM_work')
import mcore as M, engine as E
LEAVES=set(M.bools())
sd={int(k):int(v) for k,v in json.load(open('/home/user/integer_solver/solve_lab/agentM_work/deliv_seed.json')).items()}
handles={k:v for k,v in sd.items() if k not in LEAVES}
onleaves=[k for k in sd if k in LEAVES]
print('handles',len(handles),'onleaves',onleaves,flush=True)
ch=json.load(open('/home/user/integer_solver/solve_lab/agentE_work/chan_cfg0.json'))['chan']
A=ch[0]
Bside=sorted(LEAVES-set(A))
print('A',len(A),'B',len(Bside),flush=True)
def sc(seed):
    v=E.forward(seed); av=E.badatoms(v); return len(E.eqfails(av)),len(av)
base=sc(sd); print('deliverable seed in E model: fails=%d score=%d atoms=%d'%(base[0],39033-base[0],base[1]),flush=True)
random.seed(3)
As=random.sample(A,30)+[24601]
Bs=Bside
res={}
best=(10**9,None)
t0=time.time()
for i,a in enumerate(As):
    for b in Bs:
        s=dict(handles); s[a]=1; s[b]=1
        f,na=sc(s); res[(a,b)]=(f,na)
        if f<best[0]:
            best=(f,(a,b)); print('  NEW BEST fails=%d score=%d pair=(%d,%d) atoms=%d'%(f,39033-f,a,b,na),flush=True)
    print('%d/%d a=%d  best=%d %s  t=%.0fs'%(i+1,len(As),a,39033-best[0],best[1],time.time()-t0),flush=True)
pickle.dump(res,open('/home/user/integer_solver/solve_lab/agentM_work/pairscan.pkl','wb'))
print('BEST',best)
