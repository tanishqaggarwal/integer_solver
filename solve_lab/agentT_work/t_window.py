#!/usr/bin/env python3
"""AUDIT T4: reconstruct agent A's equation-closure windows A_L / R_L / K_L from agent F's
INDEPENDENT parse (39,033 atoms, not A's 42,267) and measure exactly what K_L excludes --
in particular whether ANY of the 256 selector booleans is ever a knob."""
import os,sys,pickle,json,collections
T=os.path.dirname(os.path.abspath(__file__))
F=os.path.abspath(os.path.join(T,'..','agentF_work')); sys.path.insert(0,F)
from circ2 import vars_of
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; eqrows=d['eqrows']; names=list(atoms); idx={a:i for i,a in enumerate(names)}
NA=len(names); NE=len(eqrows)
avars=[frozenset(vars_of(atoms[a])) for a in names]
v2a=collections.defaultdict(set)
for i,vs in enumerate(avars):
    for u in vs: v2a[u].add(i)
eq_atoms=[frozenset(idx[a] for k,a in row) for row in eqrows]
a2e=collections.defaultdict(set)
for e,s in enumerate(eq_atoms):
    for a in s: a2e[a].add(e)
FAIL=[12231,12270,12350,14584,18673,22044,29125]
SEL=set()
try:
    B=json.load(open(os.path.join(T,'..','agentH_work','bits.json'))); SEL=set(B['A'])|set(B['B'])
except Exception as ex: print('no bits.json',ex)
print('selector booleans known:',len(SEL))
A=set()
for e in FAIL: A|=eq_atoms[e]
print('%4s %8s %8s %8s %8s %8s %10s'%('L','atoms','eqs','vars','knobs','excluded','selectors_in_K'))
for L in range(0,25):
    R=set()
    for a in A: R|=a2e[a]
    V=set()
    for a in A: V|=avars[a]
    K=set(u for u in V if v2a[u]<=A)
    if L in (0,1,2,3,4,6,10,16,20,24):
        print('%4d %8d %8d %8d %8d %8d %10d'%(L,len(A),len(R),len(V),len(K),len(V)-len(K),len(K&SEL)))
    if L==24: break
    A2=set(A)
    for e in R: A2|=eq_atoms[e]
    A=A2
