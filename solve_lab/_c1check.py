import heal_harness as H, json, pickle
p=H.p
C=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/atomcache.pkl','rb'))
atoms=C['atoms']
def av(ai):
    a=atoms[ai]; s=0
    for vl,c in a['poly']:
        t=c
        for v in vl: t*=H.val[v]
        s+=t
    return s
d0=H.loadd('best/new_instance_partial_39013.json')
for v in range(H.NVARS): H.val[v]=d0.get(v,0)
H.forward()
print("atom1465 raw:", atoms[1465]['poly'])
print("atom29373 raw:", atoms[29373]['poly'])
print("\natom1465 val =", av(1465), " (==0:",av(1465)==0,")")
print("atom29373 val =", av(29373), " (==0:",av(29373)==0,")")
print("\nx_24468=",H.val[24468])
print("x_13682=",H.val[13682], " mod p=",H.val[13682]%p)
print("x_32989=",H.val[32989], " mod p=",H.val[32989]%p, " ==p*k:",H.val[32989]%p==0)
print("x_34243=",H.val[34243], " mod p=",H.val[34243]%p, " ==p*k:",H.val[34243]%p==0)
print("x_11436=",H.val[11436])
# So x_24468 mod p:
print("x_24468 mod p=",H.val[24468]%p)
# CONST1 from poly (the [] term)
for vl,c in atoms[1465]['poly']:
    if len(vl)==0: print("atom1465 const term c=",c, " c mod p=",c%p)
