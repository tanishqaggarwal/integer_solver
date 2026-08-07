import sys, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work/mirror')
import harness as H
S=pickle.load(open('x_seed.pkl','rb')); sd=S['seed']
known={2081,24601,6418,12553,22152,33462}
for k in sorted(sd):
    if k in known: continue
    print('### x_%d bits=%d'%(k,sd[k].bit_length()))
    for i in H.occ.get(k,[]):
        a=H.atoms[i]
        print('    [%d] %s'%(i, a if len(a)<220 else a[:220]+'...'))
