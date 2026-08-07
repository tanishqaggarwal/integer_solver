"""The instance is [k]P0 = P3 on secp256k1.  Take the one cheap shot: is k of low
Hamming weight?  Meet-in-the-middle over subsets of the 256 doubling-chain points."""
import os, sys, pickle, itertools, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gpt
from gsym2 import L, ad, P
ch=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/chain.pkl','rb'))
pts=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/bitpoints.pkl','rb'))
order=ch['order']; P0=ch['P0']
Q=[None]*256; R=P0
for i in range(256): Q[i]=R; R=gpt.add(R,R)
import g46_table as T
base=T.frame([]); B1,B2,B3=base['pts']
TARGET=B3
i2081=order.index(2081); i24601=order.index(24601)
print('chain index of bit x2081 = %d, x24601 = %d'%(i2081,i24601))
print('P1 == [2^%d]P0 ? %s ; P2 == [2^%d]P0 ? %s'%(i24601,B1==Q[i24601],i2081,B2==Q[i2081]))
print('target P3 =',TARGET)
print('P3 == [2^%d + 2^%d]P0 ? %s'%(i24601,i2081,TARGET==gpt.add(B1,B2)))
t0=time.time()
found=None
# weight 1,2
one={Q[i]:(i,) for i in range(256)}
if TARGET in one: found=one[TARGET]
S2={}
for i,j in itertools.combinations(range(256),2):
    S2[gpt.add(Q[i],Q[j])]=(i,j)
print('weight<=2 table built (%d) %.0fs'%(len(S2),time.time()-t0),flush=True)
if not found and TARGET in S2: found=S2[TARGET]
# weight 3,4 via MITM against S2 and one
if not found:
    for s,ij in S2.items():
        r=gpt.sub(TARGET,s)
        if r in one: found=ij+one[r]; break
        if r in S2:
            k=S2[r]
            if not set(k)&set(ij): found=ij+k; break
print('weight<=4 done %.0fs found=%s'%(time.time()-t0,found),flush=True)
if not found:
    S3={}
    for i,j,k in itertools.combinations(range(256),3):
        S3[gpt.add(gpt.add(Q[i],Q[j]),Q[k])]=(i,j,k)
    print('weight-3 table built (%d) %.0fs'%(len(S3),time.time()-t0),flush=True)
    if TARGET in S3: found=S3[TARGET]
    if not found:
        for s,ijk in S3.items():
            r=gpt.sub(TARGET,s)
            if r in one:
                if not set(ijk)&set(one[r]): found=ijk+one[r]; break
            if r in S2:
                k2=S2[r]
                if not set(ijk)&set(k2): found=ijk+k2; break
            if r in S3:
                k3=S3[r]
                if not set(ijk)&set(k3): found=ijk+k3; break
print('RESULT: k bits =',found,' (%.0fs)'%(time.time()-t0))
if found:
    k=sum(1<<i for i in found)
    print('k =',k)
    print('verify [k]P0 == P3 ?', gpt.mul(k,P0)==TARGET)
    json_bits=[order[i] for i in found]
    print('message boolean free inputs to set:',json_bits)
    pickle.dump({'bits':json_bits,'k':k},open('/home/user/integer_solver/solve_lab/agentG_work/solution_bits.pkl','wb'))
