import os, sys, itertools, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import g46_table as T, gpt
from gsym2 import L, ad, P
base=T.frame([])
B1=base['pts'][0]; B2=base['pts'][1]; B3=base['pts'][2]
print('base P1',B1); print('base P2',B2); print('base P3',B3); print('base D',base['D'],flush=True)
P1MOV=[47,112,542]; P2MOV=[91,438,490]
sing={}
for b in P1MOV+P2MOV:
    o=T.frame([b]); sing[b]=o
    d1=gpt.sub(o['pts'][0],B1); d2=gpt.sub(o['pts'][1],B2)
    print('x%-5d dP1=%s dP2=%s'%(b,str(d1)[:46],str(d2)[:46]),flush=True)
print('\n--- same-coordinate pairs: is the MUX ADDITIVE or EXCLUSIVE? ---')
for i,j in [(47,112),(47,542),(112,542),(91,438),(91,490),(438,490)]:
    o=T.frame([i,j])
    if 'pts' not in o: print('  (%d,%d): %s'%(i,j,o.get('core'))); continue
    which=0 if i in P1MOV else 1
    Pi=sing[i]['pts'][which]; Pj=sing[j]['pts'][which]; Bb=base['pts'][which]
    Pij=o['pts'][which]
    addpred=gpt.sub(gpt.add(Pi,Pj),Bb)
    print('  (%d,%d) nres=%d inc=%s oncurve=%s | equals bit%d:%s bit%d:%s base:%s ADDITIVE:%s'%(
        i,j,o['nres'],o.get('res_lin_inconsistent'),o.get('oncurve'),
        i,Pij==Pi,j,Pij==Pj,Bb==Pij,Pij==addpred),flush=True)
print('\n--- cross pairs: F_p-affine / group-affine on D ---')
for i,j in [(47,91),(47,438),(112,91),(542,490)]:
    o=T.frame([i,j])
    if o.get('D') is None or sing[i].get('D') is None or sing[j].get('D') is None:
        print('  (%d,%d): D undefined'%(i,j)); continue
    predD=gpt.sub(gpt.add(sing[i]['D'],sing[j]['D']),base['D'])
    print('  (%d,%d): D=%s  group-affine=%s'%(i,j,str(o['D'])[:44],o['D']==predD),flush=True)
