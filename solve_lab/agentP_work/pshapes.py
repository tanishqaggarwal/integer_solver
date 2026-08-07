#!/usr/bin/env python3
import pickle,sys,json
from collections import Counter,defaultdict
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb'))
rows,AP=D['rows'],D['AP']

def shape(ap):
    """canonical shape signature: sorted list of (monomial-pattern, coeff-class)"""
    # relabel vars by first appearance
    ren={}; out=[]
    for m,c in sorted(ap.items(), key=lambda z:(len(z[0]),z[0])):
        mm=[]
        for x in m:
            if x not in ren: ren[x]=len(ren)
            mm.append(ren[x])
        cc = c if abs(c)<=2 else ('BIG' if abs(c)>10**20 else 'C')
        if not m and abs(c)>2: cc='BIG' if abs(c)>10**20 else 'C'
        out.append((tuple(sorted(mm)),cc))
    return tuple(sorted(out))

S=Counter()
byshape=defaultdict(list)
for i,ap in enumerate(AP):
    s=shape(ap); S[s]+=1; byshape[s].append(i)
print("distinct shapes:",len(S))
for s,c in S.most_common(40):
    print(f"{c:7d}  {s}")
pickle.dump({'byshape':dict(byshape)},open(W+'shapes.pkl','wb'))
