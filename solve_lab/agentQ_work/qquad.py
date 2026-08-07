#!/usr/bin/env python3
"""Q-11b: the four quadrants of ONE slot, checked numerically against the group law.
Slot atoms, read verbatim off EQUATIONS.txt:
   a = x_2779 (sel, leaf 2^0)   b = x_34715 (sel, leaf 2^164)   both boolean-pinned
   x_3565=a  x_31966=b  x_24678=1-b  x_24849=1-a
   cA = x_13201 = x_3565*x_24678 = a(1-b)
   cB = x_33391 = x_31966*x_24849 = b(1-a)
   cC = x_4639  = x_31966*x_3565  = a*b
   Xout = x_20820 = cA*x_22231 + cB*x_11321 + cC*x_22294
   Yout = x_18440 = cA*x_27051 + cB*x_37031 + cC*x_33676
   live_out = x_11830 - x_1609 = (a+b) - ab = a OR b
"""
import json,sys
sys.path.insert(0,'.')
from qgrp import add,p,cs
leaf={int(g):v for g,v in json.load(open('qleaf.json')).items()}
A=(int(leaf[2779][0]),int(leaf[2779][1])); B=(int(leaf[34715][0]),int(leaf[34715][1]))
Xa,Ya=(A[0]-cs)%p,A[1]%p; Xb,Yb=(B[0]-cs)%p,B[1]%p
S=add(A,B); u3,y3=(S[0]-cs)%p,S[1]%p
print('slot inputs: leaf 2^0 and leaf 2^164 ; chord output = their group sum')
print('%-8s %-6s %-6s %-6s  %-28s %s'%('(a,b)','cA','cB','cC','slot carries','matches'))
for a in (0,1):
    for b in (0,1):
        cA=a*(1-b); cB=b*(1-a); cC=a*b
        X=(cA*Xa+cB*Xb+cC*u3)%p; Y=(cA*Ya+cB*Yb+cC*y3)%p
        live=(a+b)-a*b
        exp={(0,0):('identity (0,0)',(0,0)),(1,0):('leaf 2^0',(Xa,Ya)),
             (0,1):('leaf 2^164',(Xb,Yb)),(1,1):('sum 2^0 + 2^164',(u3,y3))}[(a,b)]
        print('%-8s %-6d %-6d %-6d  %-28s %s   live_out=%d'%(
            '(%d,%d)'%(a,b),cA,cB,cC,exp[0],(X,Y)==exp[1],live))
