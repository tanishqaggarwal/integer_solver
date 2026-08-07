import sys, pickle, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentL_work')
import importlib.util
spec=importlib.util.spec_from_file_location('ff','/home/user/integer_solver/solve_lab/agentL_work/fastfold.py')
ff=importlib.util.module_from_spec(spec); spec.loader.exec_module(ff)
spec2=importlib.util.spec_from_file_location('ss','/home/user/integer_solver/solve_lab/agentL_work/subsearch.py')
ss=importlib.util.module_from_spec(spec2); spec2.loader.exec_module(ss)
rnd=random.Random(99); ok=0; bad=0
for _ in range(300):
    A,B=rnd.sample(ff.live,2)
    if ff.fold([A,B])==ss.fold2(A,B): ok+=1
    else: bad+=1
print('fold2 (LCA) vs full fold on 300 random pairs: match %d mismatch %d'%(ok,bad))
