import sys, time, os
t0=time.time()
sys.argv=['x']
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentT_work')
import t_close2w as C
print('import wall %.1f s'%(time.time()-t0))
M=C.M
print('live selectors:', len(M['live']))
print('dead:', len(M['dead']))
print('NV', C.NV, 'n atoms', len(C.E.res))
print('SHIFT size', len(C.SHIFT))
print('first live', sorted(M['live'])[:10])
