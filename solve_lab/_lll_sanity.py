from fractions import Fraction as Fr
def lll(B, delta=Fr(3,4)):
    B=[row[:] for row in B]; n=len(B)
    def dot(u,v): return sum(a*b for a,b in zip(u,v))
    def gram():
        Bs=[]; mu=[[Fr(0)]*n for _ in range(n)]
        for i in range(n):
            bi=[Fr(x) for x in B[i]]
            for j in range(i):
                mu[i][j]=dot([Fr(x) for x in B[i]],Bs[j])/dot(Bs[j],Bs[j])
                bi=[bi[k]-mu[i][j]*Bs[j][k] for k in range(len(bi))]
            Bs.append(bi)
        return Bs,mu
    Bs,mu=gram(); k=1
    while k<n:
        for j in range(k-1,-1,-1):
            if abs(mu[k][j])>Fr(1,2):
                q=round(mu[k][j])
                B[k]=[B[k][t]-q*B[j][t] for t in range(len(B[k]))]
                Bs,mu=gram()
        if dot(Bs[k],Bs[k])>=(delta-mu[k][k-1]**2)*dot(Bs[k-1],Bs[k-1]): k+=1
        else:
            B[k],B[k-1]=B[k-1],B[k]; Bs,mu=gram(); k=max(k-1,1)
    return B
# textbook: should reduce to short near-orthogonal vectors
B=[[1,1,1],[-1,0,2],[3,5,6]]
R=lll(B)
print("reduced:",R)
import math
for r in R: print("  norm^2:",sum(x*x for x in r))
# known: LLL of this gives rows like [0,1,-1] etc with small norms
print("determinant preserved check (|det|): original vs reduced")
def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
print(abs(det3(B)),abs(det3(R)),"(equal => valid unimodular reduction)")
