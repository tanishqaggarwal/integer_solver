"""Measure the map  (which gating bits are on) -> (x1,y1,x2,y2, A, B)."""
import json, sys, time, collections
import dlib as L
import engine2 as E
import adv3
P = L.P
COORD = {'x1': 12186, 'y1': 16742, 'x2': 14853, 'y2': 24908, 'x3': 22162, 'y3': 30213,
         's': 9192, 'A': 35389, 'B': 6671}

t = json.load(open('table.json'))
banks = json.load(open('banks.json'))
base = L.load('D_adv.json')


def run(bits_on, bits_off, sweeps=12):
    st = E.St(base)
    seeds = {}
    for b in bits_off:
        seeds[b] = 0
    for b in bits_on:
        seeds[b] = 1
    st.apply(seeds)
    adv3.sweep(st, rounds=sweeps)
    return st


if __name__ == '__main__':
    st = E.St(base)
    print('base score', st.score, st.nz())
    print('base coords:', {k: st.v[u] % P for k, u in COORD.items()})
    b1 = banks['bank1']
    b2 = banks['bank2']
    def consts(b):
        return [C % P for a, x, C in t[str(b)]]
    print('bit 24601 consts', consts(24601))
    print('bit 2081  consts', consts(2081))
    for i in b1[:6]:
        if i == 24601:
            continue
        s = run([i], [24601])
        c = consts(i)
        print(f'bank1 bit x_{i}: score={s.score} nz={len(s.nz())} x1={s.v[12186]%P} y1={s.v[16742]%P}')
        print(f'      consts {c}  match_x1={s.v[12186]%P in c} match_y1={s.v[16742]%P in c}')
