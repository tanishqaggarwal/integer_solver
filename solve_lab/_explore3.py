import heal_harness as H, json
p=H.p
for path in ['gadget_handled.json','best_agentA_39022.json','best/new_instance_partial_39013.json']:
    try:
        d=H.loadd(path)
        for v in range(H.NVARS): H.val[v]=d.get(v,0)
        H.forward()
        F=H.fails()
        print(f"{path}: {len(F)} fail -> {F}")
    except Exception as e:
        print(f"{path}: ERROR {e}")
