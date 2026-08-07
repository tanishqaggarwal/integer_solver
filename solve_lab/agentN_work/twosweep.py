"""Re-run the cascade-pin stage-A pricing from TWO-SELECTOR bases and compare to the one-selector
   table.  Does rank rise faster than deficit when a second selector (and the accumulator) is live?"""
import ev, fast, json, time
from fast import St, csup
from chain import close_trace
from two import price          # same pricing function, unchanged
from collections import defaultdict
BITS=json.load(open('bits.json')); UW=json.load(open('uw.json'))
ALL=set(BITS['A']+BITS['B']); U=UW['U']; W=UW['W']
st0=St({})
FRAMES=[('1sel  u[0]',       {U[0]:1}),
        ('2sel  same-u',     {U[0]:1,U[1]:1}),
        ('2sel  same-u far', {U[0]:1,U[60]:1}),
        ('2sel  same-w',     {W[0]:1,W[1]:1}),
        ('2sel  witness set',{24601:1,2081:1}),
        ('3sel  u+u+w',      {U[0]:1,U[1]:1,W[0]:1})]
out={}
for name,sel in FRAMES:
    t0=time.time()
    base,ok,tr,_=close_trace(st0.clone().set_free(dict(sel)),frozen=set(ALL))
    chain=[a for a,X in tr]
    rows=[]
    for P in chain:
        s,ok2,tr2,_=close_trace(st0.clone().set_free(dict(sel)),frozen=set(ALL),skip=frozenset([P]))
        r=price(s,'a%d'%P); rows.append(r)
    base_r=price(base,'base')
    best=max(rows,key=lambda z:z['score']) if rows else base_r
    mingap=min((z['deficit']-z['rank']) for z in rows) if rows else 99
    maxrank=max((z['rank']) for z in rows) if rows else 0
    wins=[z for z in rows if z['rank']>z['deficit']]
    out[name]=dict(base_score=base.score(),chain=len(chain),best=best['score'],
                   mingap=mingap,maxrank=maxrank,wins=len(wins),
                   base_deficit=base_r['deficit'],base_rank=base_r['rank'],base_knobs=base_r['knobs'])
    print('%-18s base=%d chain=%2d | base deficit=%-3d rank=%-2d knobs=%-2d | pins: best=%d maxrank=%d min(deficit-rank)=%d  rank>deficit: %d   %.0fs'%(
        name,base.score(),len(chain),base_r['deficit'],base_r['rank'],base_r['knobs'],
        best['score'],maxrank,mingap,len(wins),time.time()-t0),flush=True)
    for z in wins: print('   *** WIN a%s rank=%d deficit=%d score=%d'%(z['tag'],z['rank'],z['deficit'],z['score']),flush=True)
json.dump(out,open('twosweep.json','w'),indent=1)
print('\nWITNESS reference (2 selectors on: x_24601, x_2081): deficit=4 rank=7 failing=7 score=39026')
