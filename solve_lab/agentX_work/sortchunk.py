import sys, numpy as np
f=sys.argv[1]; k=np.fromfile(f,dtype=np.uint64); k.sort(); k.tofile(f)
