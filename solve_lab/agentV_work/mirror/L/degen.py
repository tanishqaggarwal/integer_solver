"""Are any configurations degenerate (chord denominator 0)?  Bears on the 2^256-1 count."""
import importlib.util, itertools, random, collections
spec=importlib.util.spec_from_file_location('ss','/home/user/integer_solver/solve_lab/agentT_work/mirror/L/subsearch.py')
ss=importlib.util.module_from_spec(spec); spec.loader.exec_module(ss)
p=ss.p; live=ss.live
n=0; deg=0
for A,B in itertools.combinations(live,2):
    n+=1
    if ss.fold2(A,B) is None: deg+=1; print('DEGENERATE PAIR',A,B)
print('|S|=2: %d pairs, %d degenerate'%(n,deg))
# random larger sets via the full fold
spec2=importlib.util.spec_from_file_location('ff','/home/user/integer_solver/solve_lab/agentT_work/mirror/L/fastfold.py')
ff=importlib.util.module_from_spec(spec2); spec2.loader.exec_module(ff)
rnd=random.Random(1); bad=0
for k in (3,4,8,16,64,128,256):
    b=0
    for _ in range(60 if k<200 else 1):
        S=rnd.sample(live,k)
        if ff.fold(S) is None: b+=1
    print('  |S|=%-4d random sets: %d degenerate of %d'%(k,b,60 if k<200 else 1))
