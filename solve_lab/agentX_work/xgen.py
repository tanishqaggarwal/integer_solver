import json,os,sys
HERE=os.path.dirname(os.path.abspath(__file__))
d=json.load(open(os.path.join(HERE,'xdata.json')))
def w(path,T,lad):
    with open(path,'w') as f:
        f.write('%d %d\n'%(T[0],T[1]))
        for x,y in lad: f.write('%d %d\n'%(x,y))
lad=[(int(a),int(b)) for a,b in d['ladder']]
T=(int(d['T'][0]),int(d['T'][1]))
w(os.path.join(HERE,'data_real.txt'),T,lad)
print('wrote data_real.txt')
