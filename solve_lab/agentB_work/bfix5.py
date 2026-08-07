import os,sys,json
os.environ['ORIENT']=os.environ.get('ORIENT','orient7.pkl'); sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentB_work')
import beval as E, bfix as F
P=F.P; Q=F.Q
v0=E.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
def trial(knobs,iters=6,dump=None):
    fv={v:v0[v] for v in E.free}
    val,nd,_=E.forward(fv,default=v0)
    for it in range(iters):
        if 'A4432' in knobs: fv[4432]-= (val[28730]%P)
        if 'A12553' in knobs: fv[12553]+= (val[28730]%P)
        if 'B' in knobs: fv[3349]-= (val[8731]%P)
        if 'E' in knobs: fv[9118]-= (val[9118]%P)
        if 'D7068' in knobs: fv[7068]-= ((val[7068]-val[2099])%Q)
        if 'D6418' in knobs: fv[6418]+= ((val[7068]-val[2099])%Q)
        val,nd,_=E.forward(fv,default=v0)
    val=F.set_handles(val)
    ok,fail,fvv=E.score(val)
    r=F.residues(val)
    print(knobs,'-> score',ok,'nfail',len(fail),'nd',nd,'zeroed:',[k[0] for k,vv in r.items() if vv==0], flush=True)
    if dump and ok>=39026:
        json.dump({('x_%d'%i):val[i] for i in range(38748)},open(dump%ok,'w')); print('  WROTE',dump%ok)
    return ok,val
for ks in [['A12553'],['D6418'],['B','E'],['B','E','A12553'],['B','E','D6418'],['B','E','A12553','D6418'],
           ['B','E','A12553','D7068'],['B','E','A4432','D6418']]:
    trial(ks,dump='out/B_%d.json')
