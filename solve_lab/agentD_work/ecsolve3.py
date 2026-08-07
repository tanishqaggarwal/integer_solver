"""Solve A = B = 0 (mod p) using x_22152 -> x1 and x_33462 -> y1.

Closed form with x2,y2,x3,y3 held at their pinned values:
    B = 0  =>  y1 = [ y2 (x1 - x3) - y3 (x2 - x1) ] / (x2 - x3)
    then  y2 - y1 = (x2 - x1)(y2 + y3)/(x2 - x3)
    A = 0  =>  (x2-x1)^2 [ (x1+x2+x3+CA) - (y2+y3)^2/(x2-x3)^2 ] = 0
  branch 1: x1 = x2 (and then y1 = y2)
  branch 2: x1 = (y2+y3)^2/(x2-x3)^2 - x2 - x3 - CA
"""
import json, sys, time
import dlib as L
import engine2 as E
import adv3
import hsweep
P = L.P
X1, Y1, X2, Y2, X3, Y3 = 12186, 16742, 14853, 24908, 22162, 30213
K1, K2 = 22152, 33462          # knobs that drive x1 and y1
A_, B_ = 35389, 6671
CA = 97553848499418123410591666447050222001188385549510401465815187079080512838891

st0 = E.St(L.load(sys.argv[1] if len(sys.argv) > 1 else 'D_adv.json'))
print('base', st0.score, st0.nz())


def coords(st):
    return [st.v[u] % P for u in (X1, Y1, X2, Y2, X3, Y3)]


def run(seeds, sweeps=10):
    st = st0.clone()
    st.apply(seeds)
    adv3.sweep(st, rounds=sweeps)
    return st


c0 = coords(st0)
print('base coords', c0)
# measure the linear map knob -> (x1,y1)
p1 = coords(run({K1: st0.v[K1] + 1}))
p2 = coords(run({K2: st0.v[K2] + 1}))
print('d(x1,y1)/dK1 =', (p1[0] - c0[0]) % P, (p1[1] - c0[1]) % P)
print('d(x1,y1)/dK2 =', (p2[0] - c0[0]) % P, (p2[1] - c0[1]) % P)

x1, y1, x2, y2, x3, y3 = c0
inv = lambda z: pow(z % P, P - 2, P)
d = (x2 - x3) % P
sols = []
# branch 1
sols.append(('x1=x2', x2 % P, y2 % P))
# branch 2
if d:
    t = (y2 + y3) % P * inv(d) % P
    nx1 = (t * t - x2 - x3 - CA) % P
    ny1 = (y2 * (nx1 - x3) - y3 * (x2 - nx1)) % P * inv(d) % P
    sols.append(('quad', nx1, ny1))
for nm, nx1, ny1 in sols:
    A = ((nx1 + x2 + x3 + CA) * pow((x2 - nx1) % P, 2, P) - pow((y2 - ny1) % P, 2, P)) % P
    B = ((y3 + ny1) * (x2 - nx1) - (y2 - ny1) * (nx1 - x3)) % P
    print(f'branch {nm}: x1={nx1} y1={ny1}  A={A} B={B}')

best = None
for nm, nx1, ny1 in sols:
    st = run({K1: nx1, K2: ny1}, sweeps=12)
    print(f'  [{nm}] after sweep: score={st.score} nz={st.nz()}  A={st.v[A_]%P} B={st.v[B_]%P}')
    hsweep.sweep(st, rounds=6)
    print(f'  [{nm}] after handles: score={st.score} nz={st.nz()}')
    if best is None or st.score > best[0]:
        best = (st.score, nm, st)
print('BEST', best[0], best[1])
L.save(best[2].v, f'D_ec3_{best[0]}.json')
print('saved D_ec3_%d.json' % best[0])
