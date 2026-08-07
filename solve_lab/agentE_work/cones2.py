import sys, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import harness as H
definer=H.definer; avars=H.avars; atoms=H.atoms
def cone(aid):
    seen=set(); stack=list(avars[aid]); fr=set()
    while stack:
        u=stack.pop()
        if u in seen: continue
        seen.add(u)
        dv=definer[u]
        if dv is None: fr.add(u); continue
        for w in avars[dv[0]]:
            if w!=u: stack.append(w)
    return seen,fr
if __name__=='__main__':
    allf=set(); 
    for aid in [int(x) for x in sys.argv[1:]]:
        s,f=cone(aid)
        print(aid, atoms[aid][:70],"cone",len(s),"free",len(f))
        allf|=f
    print("union free",len(allf))
    pickle.dump(sorted(allf), open('conefree.pkl','wb'))
