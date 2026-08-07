#!/usr/bin/env python3
"""Independent verification of the restored shared tables (read-only, stdlib only)."""
import os, struct
X='/home/user/integer_solver/solve_lab/agentX_work'
t=os.path.join(X,'tbl4s.bin'); b=os.path.join(X,'bm4.bin')
st=os.path.getsize(t)
ok=[]
def chk(label,got,want):
    good = got==want
    ok.append(good)
    print('  %-22s got=%s\n  %-22s want=%s   -> %s'%(label,got,'',want,'MATCH' if good else 'MISMATCH'))
print('tbl4s.bin size = %d bytes ; size %% 8 = %d'%(st,st%8))
chk('key count', st//8, 177589056)
f=open(t,'rb')
first2=list(struct.unpack('<2Q', f.read(16)))
chk('first two keys', first2, [208528404822,231390034609])
f.seek(st-16)
last2=list(struct.unpack('<2Q', f.read(16)))
chk('last two keys', last2, [18446743699321287810,18446743880247473500])
# structural check X did not claim: strict ascending order on a stride sample
n=st//8; NS=200000; step=max(1,n//NS); prev=-1; mono=True; nsamp=0
for i in range(0,n,step):
    f.seek(i*8); v=struct.unpack('<Q',f.read(8))[0]; nsamp+=1
    if v<prev: mono=False; print('  ORDER BREAK at index %d: %d < %d'%(i,v,prev)); break
    prev=v
f.close()
print('  sorted ascending on stride sample (%d of %d keys, step %d): %s'%(nsamp,n,step,mono))
ok.append(mono)
sb=os.path.getsize(b)
print('bm4.bin size = %d bytes (expect 1<<29 = %d)'%(sb,1<<29))
ok.append(sb==(1<<29)); print('  -> %s'%('MATCH' if sb==(1<<29) else 'MISMATCH'))
print()
print('ALL CHECKS PASS' if all(ok) else 'SOME CHECKS FAILED')
raise SystemExit(0 if all(ok) else 3)
