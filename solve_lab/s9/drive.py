"""Driver: apply seeds, ripple+repair, report."""
import pickle, sys, json, collections
import harness as H
exec(open('repair.py').read().split('if __name__')[0])

def run(seeds, label, rounds=30, verbose=True, base='../best/new_instance_partial_39022.json'):
    v = H.load_assignment(base)
    ripple(v, seeds)
    ok, hist = repair_loop(v, rounds=rounds, verbose=verbose)
    nz = nonzero_atoms(v)
    print(f'[{label}] converged={ok} hist={hist} final nonzero atoms={len(nz)}: {nz[:12]}')
    return v, nz

if __name__ == '__main__':
    v0 = H.load_assignment('../best/new_instance_partial_39022.json')
    strat = sys.argv[1] if len(sys.argv)>1 else 'B'
    if strat == 'B':
        seeds = {6418: v0[7068] - 7376877*v0[642], 12553: v0[4432] - v0[28730]}
    elif strat == 'A':
        seeds = {7068: v0[2099] + 7376877*v0[642], 4432: v0[19964] + v0[28730]}
    elif strat == 'B1':   # only C1 side via x_6418
        seeds = {6418: v0[7068] - 7376877*v0[642]}
    elif strat == 'B2':   # only C2 side via x_12553
        seeds = {12553: v0[4432] - v0[28730]}
    elif strat == 'A1':
        seeds = {7068: v0[2099] + 7376877*v0[642]}
    elif strat == 'A2':
        seeds = {4432: v0[19964] + v0[28730]}
    v, nz = run(seeds, strat)
    codes,_ = H.load_equations(); fails = H.evaluate(codes, v)
    print(f'[{strat}] EQUATIONS {len(codes)-len(fails)}/{len(codes)}')
    H.save_assignment(v, f'out_{strat}.json')
