import re,glob,os
from math import comb
w=[comb(255-i,4) for i in range(256)]
tot=comb(256,5)
done=set()
for f in sorted(glob.glob('part*.log')):
    for m in re.finditer(r'i0=(\d+) done',open(f).read()): done.add(int(m.group(1)))
d=sum(w[i] for i in done)
print('i0 completed: %d/256   work done %.3e / %.3e = %.2f%%'%(len(done),d,tot,100*d/tot))
print('completed i0 set (sorted):',sorted(done)[:40],'...' if len(done)>40 else '')
