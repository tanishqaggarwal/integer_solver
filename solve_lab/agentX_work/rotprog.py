#!/usr/bin/env python3
"""Honest progress for the rotational sweep.  A partial run does NOT exhaust a fraction of
|S|=10; it excludes exactly {S : |S|=10 and S has a balanced split at some COMPLETED rotation}.
That family's measure is estimated by sampling."""
import os,random,re
done=set()
if os.path.exists('rotdone.txt'):
    for l in open('rotdone.txt'):
        m=re.match(r'ROT (\d+) DONE',l)
        if m: done.add(int(m.group(1)))
print('rotations COMPLETE (all 6 ranges DONE): %d / 128'%len(done))
if done: print('  completed:',sorted(done))
def f(S,j): return sum(1 for i in S if ((i-j)%256)<128)
random.seed(19); n=20000; cov=0
for _ in range(n):
    S=random.sample(range(256),10)
    if any(f(S,j)==5 for j in done): cov+=1
print('fraction of random |S|=10 sets already EXCLUDED: %.4f  (%d/%d sampled)'%(cov/n,cov,n))
print()
print('CORRECT partial statement:')
print('  "every |S| = 10 ON-set having a balanced 5|5 split at one of %d completed rotations is'%len(done))
print('   excluded; |S| = 10 is exhausted only when all 128 complete."')
if len(done)<128:
    print('  NOT claimable yet: "|S| <= 10 exhausted".')
