"""WIRE-CONSISTENT selector settings.

`psel.state_for` pins every non-selector free input at the witness, so turning a leaf on breaks
its two leaf pins by construction and the region explodes to ~200 rows.  That is one honest
regime, but it is not the regime the deliverable sits in.

Every leaf pin has the exact form   sel * (w - C) - m*z ,   with w a FREE input (512 of 512) and
z = a*b a product with a FREE factor (512 of 512).  So a selector setting can be made pin-exact
with no search at all:

    live leaf   ->  set its wire w := C          (so sel*(w-C) = 0)
    every leaf  ->  set a free factor of z := 0  (so m*z = 0)

That satisfies all 512 pins simultaneously for ANY live set — the pins never obstruct a selector
setting, which is U's pin result reached from my own parse.  These configurations therefore have
small regions and are directly comparable with the deliverable's.
"""
import os, sys, json, re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
import model, ev, optN
from optN import fr, FREE, BASEFV
from frameB import State
import psel

d = model.get()
atom_src = d['atom_src']
PINS = {int(k): v for k, v in json.load(open(os.path.join(HERE, 'runs', 'pinparse.json'))).items()}
PR = re.compile(r'^x_(\d+) - x_(\d+) \* x_(\d+)$')

# free factor to drive each z to zero
ZFAC = {}
for s, ps in PINS.items():
    for a, w, C, m, z in ps:
        mm = PR.match(atom_src[ev.F['definer'][z]])
        assert mm and int(mm.group(1)) == z
        f1, f2 = int(mm.group(2)), int(mm.group(3))
        ZFAC[z] = f1 if f1 in FREE else f2
        assert ZFAC[z] in FREE


def state_for(on, zero_z=True):
    fv = dict(BASEFV)
    ons = set(on)
    for s in psel.SEL:
        if s in ons:
            fv[s] = 1
        else:
            fv.pop(s, None)
    for s, ps in PINS.items():
        for a, w, C, m, z in ps:
            if s in ons:
                fv[w] = C
            if zero_z:
                fv.pop(ZFAC[z], None)
    return State(fr, fv)


def check_pins(st):
    bad = []
    for s, ps in PINS.items():
        for a, w, C, m, z in ps:
            if st.av.get(a):
                bad.append((s, a))
    return bad


if __name__ == '__main__':
    for tag, on in psel.configs():
        st = state_for(on)
        NZ, R = psel.region_of(st)
        bad = check_pins(st)
        print('%-26s live=%-4d score=%-6d |nz|=%-4d |R|=%-5d badpins=%d'
              % (tag, len(on), st.score(), len(NZ), len(R), len(bad)), flush=True)
