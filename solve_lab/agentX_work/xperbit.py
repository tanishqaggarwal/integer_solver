#!/usr/bin/env python3
"""Agent X: is there ANY per-bit information about the ON-set, from the instance or the target?
Measures the 256 leaf selectors' footprint in EQUATIONS.txt and looks for asymmetry."""
import json,os,re,collections
HERE=os.path.dirname(os.path.abspath(__file__))
d=json.load(open(os.path.join(HERE,'xdata.json')))
e2s={int(k):int(v) for k,v in d['exp2sel'].items()}
sels={v:k for k,v in e2s.items()}          # wire -> exponent
EQ=os.path.join(HERE,'..','..','EQUATIONS.txt')
lines=open(EQ).read().split('\n')
lines=[l for l in lines if l.strip()]
print('equations:',len(lines))
occ=collections.Counter(); eqs=collections.defaultdict(set)
bool_pat=collections.Counter(); unit=collections.Counter()
pat=re.compile(r'x_(\d+)')
for li,l in enumerate(lines):
    hit=set()
    for m in pat.finditer(l):
        w=int(m.group(1))
        if w in sels:
            occ[w]+=1; hit.add(w)
    for w in hit: eqs[w].add(li)
for w in sels:
    b=0
    # booleanity gadget  (x_w)*(x_w)-(x_w)  or (x_w)*((x_w)-(1))
    pass
prof=collections.Counter()
for w in sels: prof[(occ[w],len(eqs[w]))]+=1
print('distinct (occurrences, #equations) profiles over the 256 selectors:')
for k,v in sorted(prof.items()): print('   occ=%-4d eqs=%-4d : %d selectors'%(k[0],k[1],v))
# how many selectors appear in a *unit-like* equation (an equation mentioning only that selector)
solo=[w for w in sels if any(set(int(m.group(1)) for m in pat.finditer(lines[i]))=={w} for i in eqs[w])]
print('selectors appearing in an equation that mentions no other wire (=> pinned):',len(solo))
# booleanity: does every selector carry the same boolean gadget text?
btxt=collections.Counter()
for w in sels:
    s='(x_%d)*(x_%d)-(x_%d)'%(w,w,w); t='(x_%d)*(((x_%d)-(1)))'%(w,w)
    c=sum(l.count(s) for l in lines); c2=sum(l.count(t) for l in lines)
    btxt[(c>0,c2>0)]+=1
print('booleanity-gadget presence over 256 selectors:',dict(btxt))
json.dump({'occ':{str(k):v for k,v in occ.items()},'neq':{str(k):len(v) for k,v in eqs.items()}},
          open(os.path.join(HERE,'xperbit.json'),'w'))
