import os,sys,json
os.environ['ORIENT']=os.environ.get('ORIENT','orient4.pkl'); sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentB_work')
import beval as E, bfix as F
base='/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'
v0=E.load(base); P=F.P; Q=F.Q
def trial(knobs,iters=5):
    freeval={v:v0[v] for v in E.free}
    val,nd,_=E.forward(freeval,default=v0)
    for it in range(iters):
        if 4432 in knobs: freeval[4432]=freeval[4432]-(val[28730]%P)
        if 3349 in knobs: freeval[3349]=freeval[3349]-(val[8731]%P)
        if 9118 in knobs: freeval[9118]=freeval[9118]-(val[9118]%P)
        if 7068 in knobs: freeval[7068]=freeval[7068]-((val[7068]-val[2099])%Q)
        val,nd,_=E.forward(freeval,default=v0)
    val=F.set_handles(val)
    ok,fail,fv=E.score(val)
    return ok,fail,F.residues(val),nd,val
for ks in [[],[4432],[3349],[9118],[7068],[3349,9118],[4432,3349,9118],[4432,3349,9118,7068]]:
    ok,fail,r,nd,val=trial(ks)
    print(ks,'-> score',ok,'nfail',len(fail),'nd',nd,'zeroed:',[k[0] for k,v in r.items() if v==0], flush=True)
    if ok>=39026:
        p='out/cand_%d_%s.json'%(ok,'_'.join(map(str,ks)))
        json.dump({('x_%d'%i):val[i] for i in range(38748)},open(p,'w')); print('  wrote',p)
