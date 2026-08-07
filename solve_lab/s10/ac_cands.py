"""S12 step 10: choose the TARGETED activation sets for the closure rebuild."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad, ac_lib as A
D = json.load(open(os.path.join(HERE,'ac_single.json')))
pool = D['pool']; S = {int(z): r for z, r in D['res'].items()}
zc  = [z for z in pool if S[z][0]['lost'] == 0]                       # cost 0 eqs
knob= sorted((S[z][0]['knobs'], -S[z][0]['lost'], z) for z in pool if S[z][0]['knobs'] > 0)
knob.sort(key=lambda t: (-t[0], -t[1]))
cheapknob = [t[2] for t in knob][:8]
print('zero-cost singles', len(zc), zc)
print('best knob singles (knobs, -lost, x):', knob[:10])
# free inputs that break NO check atom (they may still add closure columns)
nochk = [z for z in pool if len(S[z][0]['newchk']) == 0]
print('singles breaking no check atom:', len(nochk))
cands = [{'tag':'base','act':{}}]
for z in cheapknob[:4]:
    cands.append({'tag': f'single_{z}', 'act': {str(z): '1'}})
cands.append({'tag':'zerocost_all', 'act': {str(z): '1' for z in zc}})
cands.append({'tag':'nochk_all', 'act': {str(z): '1' for z in nochk}})
cands.append({'tag':'zerocost+24365', 'act': dict({str(z):'1' for z in zc}, **{'24365':'1'})})
cands.append({'tag':'top4knob', 'act': {str(z):'1' for z in cheapknob[:4]}})
cands.append({'tag':'allknob', 'act': {str(t[2]):'1' for t in knob}})
cands.append({'tag':'everything', 'act': {str(z):'1' for z in pool}})
json.dump(cands, open(os.path.join(HERE,'ac_cands.json'),'w'))
print('wrote', len(cands), 'candidate activation sets')
