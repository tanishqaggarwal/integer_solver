#!/usr/bin/env python3
"""AUDIT T25b -- count S's atom terms and cross-reference against F's certified parse."""
import re,pickle,sys,os,collections
LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0,os.path.join(LAB,'agentF_work'))
lhs=open('/home/user/integer_solver/EQUATIONS.txt').read().split('\n')[8680].rsplit('=',1)[0].strip()
def strip_outer(s):
    while s.startswith('(') and s.endswith(')'):
        d=0; ok=True
        for k,ch in enumerate(s):
            if ch=='(': d+=1
            elif ch==')':
                d-=1
                if d==0 and k!=len(s)-1: ok=False; break
        if not ok: break
        s=s[1:-1]
    return s
def split_top(s,op):
    s=strip_outer(s); out=[]; d=0; cur=''
    for ch in s:
        if ch=='(': d+=1
        elif ch==')': d-=1
        if d==0 and ch==op: out.append(cur); cur=''; continue
        cur+=ch
    out.append(cur); return out
S=strip_outer(split_top(strip_outer(split_top(lhs,'*')[0]),'*')[0])
# flatten the left-nested + chain
def flat(s):
    parts=split_top(s,'+')
    if len(parts)==1: return [s]
    return flat(parts[0])+parts[1:]
terms=flat(S)
print('S flattens to %d top-level + terms'%len(terms))
coefs=[]
for t in terms:
    t=t.strip()
    m=re.match(r'^\((-?\d+)\)\*\((.*)\)$',t)
    if m: coefs.append((int(m.group(1)),strip_outer(m.group(2))))
    else: coefs.append((1,strip_outer(t)))
print('   coefficients: %s'%[c for c,_ in coefs])
print('\n   the %d atom terms of S:'%len(coefs))
for c,a in coefs:
    print('     %+4d  %s'%(c,re.sub(r'\s','',a)[:70]))
# cross-check against F's certified-faithful parse
d=pickle.load(open(os.path.join(LAB,'agentF_work','circ4.pkl'),'rb'))
atoms=d['atoms']; eqrows=d['eqrows']
Frow=eqrows[8680]
print('\nF\'s certified parse of equation 8680: %d (coef,atom) entries'%len(Frow))
print('   F coefficients: %s'%sorted(k for k,_ in Frow))
Fatoms={re.sub(r'\s','',a) for _,a in Frow}
Ours={re.sub(r'[\s()]','',a) for _,a in coefs}
Fn={re.sub(r'[\s()]','',a) for a in Fatoms}
print('   F atom set == my S atom set (paren-insensitive)? %s'%(Fn==Ours))
print('   in F not in S: %s'%sorted(Fn-Ours)[:3])
print('   in S not in F: %s'%sorted(Ours-Fn)[:3])
print('\n== cross-link to my sixth-pass finding ==')
for u in ('x_10422','x_15120','x_35531'):
    print('   %s appears in eq8680: %s'%(u,u in lhs))
