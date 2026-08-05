#!/usr/bin/env python3
import heal_harness as H
# Do the loads / core quantities depend on x_4432 or x_7068?
core_targets={'S=x_35389':35389,'T=x_6671':6671,'L1=x_11150':11150,'L2=x_25739':25739,
              'L3=x_37758':37758,'x_15298':15298,'x_29804':29804,'x_35605':35605,
              'x_4007':4007,'x_24908':24908,'x_29322':29322,'x_3558':3558,'x_27713':27713}
for lbl,t in core_targets.items():
    a=H.anc.get(t,set())
    has=('x_4432' if 4432 in a else '')+(' x_7068' if 7068 in a else '')
    print(f"{lbl}: depends on 4432/7068? [{has.strip() or 'NO'}]  (#free anc={len(a)})")
# Reverse: what do x_4432, x_7068 feed? (descendant gate outputs) and are any core?
for w in [4432,7068]:
    desc=[t for t in H.order if w in H.anc.get(t,set())]
    core_desc=[t for t in desc if t in set(core_targets.values())]
    print(f"\nx_{w} feeds {len(desc)} gate outputs; core among them: {core_desc}")
    print(f"   sample descendants: {desc[:15]}")
