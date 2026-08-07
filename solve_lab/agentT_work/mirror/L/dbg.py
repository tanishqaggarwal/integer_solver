import sys
src=open('mkassign2.py').read().split('#MAINSTART')[0]
exec(src)
S=[M['live'][0]]
v,isl,valn=assignment(set(S),ORIENT); v[24468]=T1; v[18956]=T2
vv=[0]*NV
for k,x in v.items(): vv[k]=x
bad=relift(vv); print('round0 bad',len(bad))
r=E.run(vv)
for a in bad:
    i=E.residx[a]; cur=r[i]; sm=abs(SL[a])
    print('ATOM',a[:80],'| cur%p==0:',cur%p==0,'| cur%sm==0:',cur%sm==0)
    cands=[x for x in atomvalvars[a] if x in SHIFT]
    print('   shiftable cands',cands)
    if cur%p: print('   skipped: not 0 mod p'); continue
    for w in cands:
        old=vv[w]; vv[w]=old+p; d=E.run(vv)[i]-cur; vv[w]=old
        print('   w=%d  d/p=%s'%(w,d//p if d%p==0 else d))

print('--- apply fix manually ---')
a=[x for x in bad if 'x2055' in x][0]
i=E.residx[a]; cur=E.run(vv)[i]; sm=abs(SL[a]); w=34825
old=vv[w]; vv[w]=old+p; d=E.run(vv)[i]-cur; vv[w]=old
g=gcd(d,sm); mm=sm//g
t=(-(cur//g))*pow((d//g)%mm,-1,mm)%mm if mm>1 else 0
print('sm/p',sm//p,'d/p',d//p,'g/p',g//p,'mm',mm,'t',t)
vv[w]=old+p*t
new=E.run(vv)[i]; print('new atom value % sm ==0:',new%sm==0)
b2=relift(vv); print('unliftable after:',len(b2))
r=E.run(vv); print('nonzero atoms:',sum(1 for x in r if x))
for x in b2: print('   still bad:',x[:140])
