import numpy as np, os, time
t0=time.time()
k=np.fromfile('tbl4.bin',dtype=np.uint64)
print('loaded',len(k),'%.1fs'%(time.time()-t0))
assert len(k)==177589056
k.sort(kind='stable')
print('sorted %.1fs'%(time.time()-t0))
k.tofile('tbl4s.bin')
print('done')
