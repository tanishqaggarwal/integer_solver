import sys, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work/mirror')
import engine3 as E3, harness as H
eng=E3.Eng(E3.BASE_DEMOTE)
FREE=set(eng.FREE)
S=pickle.load(open('x_seed.pkl','rb')); sd=S['seed']
RX=[2498,2964,6083,7068,12186,14853,22152,22649,27982]
def show(u,tag=''):
    print('--- x_%d %s free=%s def=%s'%(u,tag,u in FREE, eng.definer[u] is not None and eng.definer[u][0]))
    for i in H.occ.get(u,[]):
        a=H.atoms[i]
        if len(a)>150: a=a[:150]+'...'
        print('   [%d]%s %s'%(i,'*' if E3.ATOM2VAR.get(i)==u else ' ',a))
for u in RX: show(u,'(X-route)')
