import sys,os,json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from qgrp import *
lad=json.load(open('ladder.json'))
pins=json.load(open(os.path.join(HERE,'..','agentF_work','pins.json')))
short={g:v for g,v in pins.items() if len(v)!=2}
print('pins with only one constant:', list(short))
G=leaves()[int(lad['G_leafvar'])]
present={int(i):v for i,v in lad['ladder'].items()}
for i in lad['missing']:
    Q=mul(pow(2,i),G); raw=((Q[0]-cs)%p, Q[1]%p)
    print('exponent %d predicted raw (x,y) = %s'%(i,raw))
    for g,v in short.items():
        for w,c in v:
            if c%p==raw[0]%p: print('    MATCH x  <- pin selector %s wire %d'%(g,w))
            if c%p==raw[1]%p: print('    MATCH y  <- pin selector %s wire %d'%(g,w))
