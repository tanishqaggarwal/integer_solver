#!/usr/bin/env python3
"""Read-only: check recorded scan counts against closed forms, and check AA's
shell-written SHARD markers against engine-written DONE lines."""
import os,re,glob
from math import comb
L='/home/user/integer_solver/solve_lab'

print('=== agent Y: rep_comp.txt counts vs C(256,b) ===')
p=os.path.join(L,'agentY_work','rep_comp.txt')
for line in open(p):
    m=re.match(r'DONE size=(\d+) range=\[(\d+),(\d+)\) n=(\d+)',line)
    if m:
        b,lo,hi,n=(int(x) for x in m.groups())
        c=comb(hi-lo,b)
        print('  size=%d range=[%d,%d) n=%-12d C(%d,%d)=%-12d %s'
              %(b,lo,hi,n,hi-lo,b,c,'MATCH' if n==c else '*** MISMATCH ***'))

print()
print('=== agent AA: shell SHARD markers vs engine DONE lines, per runs file ===')
bad=[]
for f in sorted(glob.glob(os.path.join(L,'agentAA_work','runs*','*.txt'))):
    txt=open(f,errors='replace').read().splitlines()
    sh=sum(1 for l in txt if re.match(r'^SHARD\d+$',l))
    dn=sum(1 for l in txt if l.startswith('DONE '))
    if sh==0 and dn==0: continue
    flag=''
    if sh>dn:
        flag='  *** %d MARKER(S) WITH NO ENGINE DONE ***'%(sh-dn); bad.append(f)
    print('  %-46s SHARD=%-3d DONE=%-3d%s'%(os.path.basename(f),sh,dn,flag))
print()
print('files where a shell marker outran the engine: %d'%len(bad))
for f in bad: print('   ',f)

print()
print('=== agent AA: within a tag+size, all 8 table-shards must report the same n ===')
for f in sorted(glob.glob(os.path.join(L,'agentAA_work','runs','rs_*.txt'))):
    groups={}
    for l in open(f,errors='replace'):
        m=re.match(r'DONE \S+ sz=(\d+) range=\[(\d+),(\d+)\) n=(\d+)',l)
        if m:
            b,lo,hi,n=(int(x) for x in m.groups())
            groups.setdefault((b,lo,hi),[]).append(n)
    out=[]
    for k in sorted(groups):
        v=groups[k]
        out.append('sz=%d:%s(x%d)'%(k[0],'CONSISTENT' if len(set(v))==1 else '*** SPLIT %s ***'%sorted(set(v)),len(v)))
    if out: print('  %-28s %s'%(os.path.basename(f),' '.join(out)))
