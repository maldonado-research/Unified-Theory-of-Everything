#!/usr/bin/env python3
"""Exact curvature/derivative tensor ranks and sourced placement audit.

Standard-library Gaussian rational arithmetic. No repository imports,
floating-point arithmetic, Python assertions, or loop-residue inputs.
"""
from collections import Counter
from fractions import Fraction as R
from itertools import combinations, product, permutations
import json
from pathlib import Path
import random
import sys


class Q:
    __slots__=('a','b')
    def __init__(self,a=0,b=0):
        if isinstance(a,Q):self.a,self.b=a.a,a.b
        else:self.a,self.b=R(a),R(b)
    def __add__(self,o):o=Q(o);return Q(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return Q(-self.a,-self.b)
    def __sub__(self,o):return self+-Q(o)
    def __rsub__(self,o):return Q(o)+-self
    def __mul__(self,o):
        if isinstance(o,Grass):return o*self
        o=Q(o);return Q(self.a*o.a-self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=Q(o);d=o.a*o.a+o.b*o.b
        return Q((self.a*o.a+self.b*o.b)/d,(self.b*o.a-self.a*o.b)/d)
    def __eq__(self,o):o=Q(o);return self.a==o.a and self.b==o.b
    def __bool__(self):return bool(self.a or self.b)
    def conj(self):return Q(self.a,-self.b)
    def __repr__(self):return 'Q(%s,%s)'%(self.a,self.b)


class Grass:
    def __init__(self,d=None):self.d={k:Q(v) for k,v in (d or {}).items() if v}
    @staticmethod
    def cast(x):return x if isinstance(x,Grass) else Grass({():Q(x)})
    def __add__(self,o):
        o=Grass.cast(o);d=dict(self.d)
        for k,v in o.d.items():d[k]=d.get(k,Q())+v
        return Grass(d)
    __radd__=__add__
    def __neg__(self):return Grass({k:-v for k,v in self.d.items()})
    def __sub__(self,o):return self+-Grass.cast(o)
    def __mul__(self,o):
        o=Grass.cast(o);d={}
        for a,x in self.d.items():
            for b,y in o.d.items():
                if set(a)&set(b):continue
                key=tuple(sorted(a+b));sign=(-1)**sum(i>j for i in a for j in b)
                d[key]=d.get(key,Q())+sign*x*y
        return Grass(d)
    __rmul__=__mul__
    def dagger(self):
        out=Grass()
        for k,v in self.d.items():
            term=Grass.cast(v.conj())
            for j in reversed(k):term=term*Grass({(j^1,):1})
            out+=term
        return out
    def __eq__(self,o):return self.d==Grass.cast(o).d
    def __repr__(self):return repr(self.d)


I=Q(0,1);eta=(1,-1,-1,-1)
def matrix(n,m,f):return tuple(tuple(Q(f(i,j)) for j in range(m)) for i in range(n))
def eye(n):return matrix(n,n,lambda i,j:i==j)
def zero(n,m=None):return matrix(n,n if m is None else m,lambda i,j:0)
def scale(c,a):return tuple(tuple(c*x for x in row) for row in a)
def add(a,b):return tuple(tuple(x+y for x,y in zip(ar,br)) for ar,br in zip(a,b))
def sub(a,b):return add(a,scale(-1,b))
def mul(a,b):return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(len(b))),Q()) for j in range(len(b[0]))) for i in range(len(a)))
def transpose(a):return tuple(zip(*a))
def dagger(a):return tuple(tuple(a[j][i].conj() for j in range(len(a))) for i in range(len(a[0])))
def comm(a,b):return sub(mul(a,b),mul(b,a))
def msum(xs,n=4,m=None):
    out=zero(n,m)
    for x in xs:out=add(out,x)
    return out
def tr(a):return sum((a[i][i] for i in range(len(a))),Q())
def eps(a,b,c):
    if len({a,b,c})<3:return 0
    return 1 if (a,b,c) in ((0,1,2),(1,2,0),(2,0,1)) else -1
def eps4(*v):return 0 if len(set(v))<4 else (-1)**sum(v[i]>v[j] for i in range(4) for j in range(i+1,4))
pauli=(matrix(2,2,lambda i,j:i!=j),((Q(),-I),(I,Q())),((Q(1),Q()),(Q(),Q(-1))))
T=tuple(scale(R(1,2),p) for p in pauli)
g0=matrix(4,4,lambda i,j:(1 if i<2 else -1) if i==j else 0)
def spatial(p):return tuple(tuple(p[i][j-2] if i<2 and j>=2 else -p[i-2][j] if i>=2 and j<2 else Q() for j in range(4)) for i in range(4))
gamma=(g0,)+tuple(spatial(p) for p in pauli)
sigma=tuple(tuple(scale(I/2,comm(gamma[m],gamma[n])) for n in range(4)) for m in range(4))
g5=scale(I,mul(mul(mul(gamma[0],gamma[1]),gamma[2]),gamma[3]))
PL=scale(R(1,2),sub(eye(4),g5));PR=scale(R(1,2),add(eye(4),g5));C=scale(I,mul(gamma[2],gamma[0]))
def grade3(m,n,r):
    labels=(m,n,r)
    return msum((scale(R((-1)**sum(p[a]>p[b] for a in range(3) for b in range(a+1,3)),6),mul(mul(gamma[labels[p[0]]],gamma[labels[p[1]]]),gamma[labels[p[2]]])) for p in permutations(range(3))))
G=tuple(tuple(tuple(grade3(m,n,r) for r in range(4)) for n in range(4)) for m in range(4))
def dual(ff):
    return [[msum((scale(R(eta[m]*eta[n]*eps4(m,n,r,s),2),ff[r][s]) for r in range(4) for s in range(4)),len(ff[0][0])) for n in range(4)] for m in range(4)]


counts=Counter();mutations=set()
def check(label,a,b):
    if a!=b:raise RuntimeError('check failed: '+label+'; '+repr(a)+' != '+repr(b))
    counts[label]+=1
def reject(label,a,b):
    if a==b:raise RuntimeError('mutation survived: '+label)
    mutations.add(label)
def rank(rows,ncols=None):
    if ncols is None:ncols=len(rows[0]) if rows else 0
    a=[list(row) for row in rows];r=0
    for c in range(ncols):
        p=next((i for i in range(r,len(a)) if a[i][c]),None)
        if p is None:continue
        a[r],a[p]=a[p],a[r];v=a[r][c];a[r]=[x/v for x in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][c]:
                z=a[i][c];a[i]=[x-z*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r
def flat(a):return tuple(x for row in a for x in row)


for m,n,r in product(range(4),repeat=3):
    expected=msum((scale(I*eta[s]*eps4(m,n,r,s),mul(gamma[s],g5)) for s in range(4)))
    check('every triple-gamma dual coefficient',G[m][n][r],expected)
check('chiral vector bilinear spans four vector components',rank([flat(mul(mul(PR,gamma[m]),PL)) for m in range(4)],16),4)


# Complete Lorentz invariant problem: a general coefficient of
# V^a w^b F^{rs}, r<s, has 96 unknowns. No tensor ansatz is supplied.
pairs=tuple(combinations(range(4),2));keys=tuple((a,b,r,s) for a,b in product(range(4),repeat=2) for r,s in pairs);idx={k:i for i,k in enumerate(keys)}
gens=[]
for r,s in pairs:
    gens.append(matrix(4,4,lambda a,b:1 if (a,b)==(r,s) else (-eta[r]*eta[s] if (a,b)==(s,r) else 0)))
rows=[[Q() for k in keys] for gen in gens for row in keys]
for c,(a,b,r,s) in enumerate(keys):
    for gi,gen in enumerate(gens):
        off=gi*len(keys)
        for v in range(4):
            rows[off+idx[v,b,r,s]][c]+=gen[a][v]
            rows[off+idx[a,v,r,s]][c]+=gen[b][v]
            if v!=s:
                lo,hi=sorted((v,s));sign=1 if v<s else -1
                rows[off+idx[a,b,lo,hi]][c]+=sign*gen[r][v]
            if r!=v:
                lo,hi=sorted((r,v));sign=1 if r<v else -1
                rows[off+idx[a,b,lo,hi]][c]+=sign*gen[s][v]
lorank=rank(rows,len(keys));check('full Lorentz covariance nullity',len(keys)-lorank,2)
ordinary=tuple(Q(eta[a]*eta[b]*(int(a==r and b==s)-int(a==s and b==r))) for a,b,r,s in keys)
dualvec=tuple(Q(-eps4(a,b,r,s)) for a,b,r,s in keys)
for v in (ordinary,dualvec):check('explicit ordinary or dual solves all Lorentz equations',tuple(sum((x*y for x,y in zip(row,v)),Q()) for row in rows),tuple(Q() for row in rows))
check('ordinary and dual span complete Lorentz nullspace',rank([ordinary,dualvec],len(keys)),2)


weak_receipts=[]
for weak in (False,True):
    wk=tuple((j,a,b) for j in range(3 if weak else 1) for a,b in product(range(2),repeat=2));wi={k:i for i,k in enumerate(wk)}
    wr=[[Q() for k in wk] for gen in range(3) for row in wk]
    for c,(j,a,b) in enumerate(wk):
        for gen,t in enumerate(T):
            off=gen*len(wk)
            for v in range(2):
                wr[off+wi[j,a,v]][c]+=t[b][v]
                wr[off+wi[j,v,b]][c]-=t[v][a]
            if weak:
                for jout in range(3):wr[off+wi[jout,a,b]][c]+=I*eps(gen,jout,j)
    v=tuple(T[j][a][b] if weak else Q(a==b) for j,a,b in wk)
    rk=rank(wr,len(wk));check('complete derivative weak intertwiner nullity',len(wk)-rk,1)
    check('explicit weak intertwiner solves every equation',tuple(sum((x*y for x,y in zip(row,v)),Q()) for row in wr),tuple(Q() for row in wr))
    weak_receipts.append(dict(strength='W' if weak else 'B',unknowns=len(wk),equations=len(wr),rank=rk,nullity=len(wk)-rk))


# Exact circular Fourier mode test. k,e are contravariant and the outgoing
# mode is exp(+i k.x). No physical external-state conjugation is assumed.
k=(Q(1),Q(),Q(),Q(1));circular=[]
for sign in (1,-1):
    e=(Q(),Q(1),sign*I,Q())
    ff=[[((I*(eta[m]*k[m]*eta[n]*e[n]-eta[n]*k[n]*eta[m]*e[m]),),) for n in range(4)] for m in range(4)]
    dd=dual(ff)
    for m,n in pairs:check('circular Fourier mode dual eigenvalue',dd[m][n],scale(-sign*I,ff[m][n]))
    plus=[[scale(R(1,2),add(ff[m][n],scale(I,dd[m][n]))) for n in range(4)] for m in range(4)]
    for m,n in pairs:check('selected selfdual Fourier channel',plus[m][n],ff[m][n] if sign==1 else zero(1))
    circular.append(dict(polarization_y='i' if sign==1 else '-i',dual_eigenvalue='-i' if sign==1 else '+i',selected_channel_nonzero=sign==1))
reject('dual sign selects opposite circular channel',circular[0]['selected_channel_nonzero'],False)


source_image_receipts=[]
for na,ns in ((1,1),(2,3)):
    p=na*ns;gp=R(3,5);g=R(2,3);yl=R(-1,2)
    blocks=[scale(w,eye(p)) for w in (gp*yl/2,I*gp*yl/2,g/2,I*g/2)]
    source_rows=[row for block in blocks for row in block]
    check('selected algebraic coefficient image has one matrix freedom',rank(source_rows,p),p)
    check('B dual coefficient is positive i ordinary',blocks[1],scale(I,blocks[0]))
    check('W dual coefficient is positive i ordinary',blocks[3],scale(I,blocks[2]))
    check('B and W selected coefficient matrices are correlated',scale(gp*yl,blocks[2]),scale(g,blocks[0]))
    check('normalized B plus channel coefficient',sub(blocks[0],scale(I,blocks[1])),scale(gp*yl,eye(p)))
    check('normalized W plus channel coefficient',sub(blocks[2],scale(I,blocks[3])),scale(g,eye(p)))
    source_image_receipts.append(dict(active_flavors=na,sterile_flavors=ns,algebraic_coordinates=4*p,raw_F_source_image_rank=p))


rng=random.Random(2026101612)
def real():return R(rng.randint(-3,3),rng.randint(1,3))
def scalar():return Q(real(),real())
def rm(n,m):return matrix(n,m,lambda a,b:scalar())
def row(P=PR):return mul(rm(1,4),P)
def col(P=PL):return mul(P,rm(4,1))
def herm():return msum((scale(real(),p) for p in (eye(2),)+T),2)
def aw(a,v):return tuple(sum((a[i][j]*v[j] for j in range(2)),Q()) for i in range(2))
def rw(l,a):return tuple(msum((scale(a[b][c],l[b]) for b in range(2)),1,4) for c in range(2))
def bil(l,ph,n,s=eye(4)):return sum((ph[a]*mul(mul(l[a],s),n)[0][0] for a in range(2)),Q())


for trial in range(8):
    gp=R(3,5);g=R(2,3);yl=R(-1,2);na,ns=2,3;f=rm(na,ns)
    aa=[herm() for m in range(4)];da=[[herm() for n in range(4)] for m in range(4)]
    dda={(m,n,r):herm() for m in range(4) for n in range(4) for r in range(n,4)}
    dd=lambda m,n,r:dda[m,min(n,r),max(n,r)]
    ff=[[add(sub(da[n][m],da[m][n]),scale(I,comm(aa[m],aa[n]))) for n in range(4)] for m in range(4)]
    df=[[[add(add(sub(dd(n,m,r),dd(m,n,r)),scale(I,add(comm(da[m][r],aa[n]),comm(aa[m],da[n][r])))),scale(I,comm(aa[r],ff[m][n]))) for n in range(4)] for m in range(4)] for r in range(4)]
    ft=dual(ff);ftt=dual(ft)
    for m,n in pairs:check('Lorentzian dual squared is minus one',ftt[m][n],scale(-1,ff[m][n]))
    phi=tuple(scalar() for a in range(2));dp=[tuple(scalar() for a in range(2)) for m in range(4)]
    dpc=[tuple(dp[m][a]+I*aw(aa[m],phi)[a] for a in range(2)) for m in range(4)]
    ls=[tuple(row() for a in range(2)) for al in range(na)]
    dl=[[tuple(row() for a in range(2)) for m in range(4)] for al in range(na)]
    dc=[[tuple(sub(dl[al][m][a],scale(I,rw(ls[al],aa[m])[a])) for a in range(2)) for m in range(4)] for al in range(na)]
    nl=[col() for i in range(ns)];nr=[col(PR) for i in range(ns)]
    dn=[[col() for m in range(4)] for i in range(ns)]
    an=[msum((scale(I,mul(gamma[m],dn[i][m])) for m in range(4)),4,1) for i in range(ns)]
    kl=[tuple(msum((scale(-I,mul(dc[al][m][a],gamma[m])) for m in range(4)),1,4) for a in range(2)) for al in range(na)]
    sb=[tuple(row(PL) for a in range(2)) for al in range(na)];sr=[col(PR) for i in range(ns)]
    mass=matrix(ns,ns,lambda i,j:(2,3,5)[i] if i==j else 0)
    rr=[add(msum((scale(mass[i][j],nr[j]) for j in range(ns)),4,1),sr[i]) for i in range(ns)]
    er=[sub(an[i],rr[i]) for i in range(ns)];eb=[tuple(sub(kl[al][a],sb[al][a]) for a in range(2)) for al in range(na)]
    def b(rows,ph,z,s=eye(4)):return sum((f[al][i]*bil(rows[al],ph,z[i],s) for al in range(na) for i in range(ns)),Q())
    def tensor(rows,z):return sum((b(rows,aw(ff[m][n],phi),z,sigma[m][n]) for m in range(4) for n in range(4)),Q())
    h=Q();ht=Q();qer=Q();qvr=Q();qeh=Q();qvh=Q();qvp=Q();qep=Q();qsigma=Q();qdiv=Q();boundary=Q();boundaryv=Q();boundarye=Q()
    for m,n in product(range(4),repeat=2):
        h+=eta[n]*b(ls,aw(ff[m][n],dpc[n]),nl,gamma[m]);ht+=eta[n]*b(ls,aw(ft[m][n],dpc[n]),nl,gamma[m])
        qvh-=eta[n]*b(ls,aw(ff[m][n],dpc[n]),nl,gamma[m])
        qvr-=eta[n]*b([dc[al][n] for al in range(na)],aw(ff[m][n],phi),nl,gamma[m])
        qdiv+=eta[n]*b(ls,aw(df[n][n][m],phi),nl,gamma[m])
        qvp+=eta[n]*b(ls,aw(ff[m][n],phi),[dn[i][n] for i in range(ns)],gamma[m])
        boundaryv+=eta[n]*(b([dc[al][n] for al in range(na)],aw(ff[m][n],phi),nl,gamma[m])+b(ls,aw(df[n][m][n],phi),nl,gamma[m])+b(ls,aw(ff[m][n],dpc[n]),nl,gamma[m])+b(ls,aw(ff[m][n],phi),[dn[i][n] for i in range(ns)],gamma[m]))
        for r in range(4):
            s=mul(sigma[m][n],gamma[r])
            drows=[dc[al][r] for al in range(na)]
            first=b(drows,aw(ff[m][n],phi),nl,s)
            qsigma+=I*first/4
            qeh-=b(ls,aw(ff[m][n],dpc[r]),nl,G[m][n][r])/2
            qer-=b(drows,aw(ff[m][n],phi),nl,G[m][n][r])/2
            qep+=b(ls,aw(ff[m][n],phi),[dn[i][r] for i in range(ns)],G[m][n][r])/2
            boundarye+=(b(drows,aw(ff[m][n],phi),nl,G[m][n][r])+b(ls,aw(df[r][m][n],phi),nl,G[m][n][r])+b(ls,aw(ff[m][n],dpc[r]),nl,G[m][n][r])+b(ls,aw(ff[m][n],phi),[dn[i][r] for i in range(ns)],G[m][n][r]))/2
            boundary+=I*(first+b(ls,aw(df[r][m][n],phi),nl,s)+b(ls,aw(ff[m][n],dpc[r]),nl,s)+b(ls,aw(ff[m][n],phi),[dn[i][r] for i in range(ns)],s))/4
    hplus=(h+I*ht)/2
    check('chiral triple gamma gives the positive i dual source',qeh,-I*ht)
    check('actual selected Higgs derivative source is one dual channel',-(qvh+qeh)/2,hplus)
    check('exact row Clifford equation before EOM',qer,qvr-tensor(kl,nl)/2)
    check('exact sterile Clifford equation before EOM',qvp+qep,-tensor(ls,an)/2)
    check('ordinary IBP retains row scalar divergence sterile placements',qvp,qvr+qvh+qdiv+boundaryv)
    check('triple-gamma IBP retains both matter derivative placements',qep,qer+qeh+boundarye)
    check('two vector boundaries reproduce compact Clifford boundary',(boundaryv+boundarye)/2,-boundary)
    check('literal OD equals row route before EOM',qsigma+qdiv/3,qvr-tensor(kl,nl)/4+qdiv/3)
    check('full placement IBP with explicit boundary',hplus,qvr-tensor(kl,nl)/4+qdiv/2+tensor(ls,an)/4-boundary)
    div=[msum((scale(eta[n]*eta[m],df[n][n][m]) for n in range(4)),2) for m in range(4)]
    jy=[real() for m in range(4)];jw=[[real() for m in range(4)] for a in range(3)]
    eB=[tr(div[m])/(2*gp*yl)-gp*jy[m] for m in range(4)]
    eW=[[2*tr(mul(T[a],div[m]))/g-g*jw[a][m] for m in range(4)] for a in range(3)]
    current=Q();gres=Q()
    for m in range(4):
        j=eta[m]*b(ls,phi,nl,gamma[m]);ji=[eta[m]*b(ls,aw(t,phi),nl,gamma[m]) for t in T]
        current+=gp*gp*yl*j*jy[m]+g*g*sum((ji[a]*jw[a][m] for a in range(3)),Q())
        gres+=gp*yl*j*eB[m]+g*sum((ji[a]*eW[a][m] for a in range(3)),Q())
    check('full matrix divergence splits into gauge currents and residuals',qdiv,current+gres)
    alt=qvr-tensor(sb,nl)/4+current/2+tensor(ls,rr)/4
    residual=-tensor(eb,nl)/4+gres/2+tensor(ls,er)/4
    check('complete sourced curvature placement certificate',hplus,alt+residual-boundary)
    # The rest of the actual P_red contains -current/6 and -tensor R/12.
    combined=hplus-current/6-tensor(ls,rr)/12
    alternate=qvr-tensor(sb,nl)/4+current/3+tensor(ls,rr)/6
    check('whole selected packet correlated coefficients',combined,alternate+residual-boundary)
    if trial==0:
        reject('opposite chiral dual sign',qeh,I*ht)
        reject('delete the transferred current',hplus,alt+residual-boundary-current/2)
        reject('delete lepton tensor source',hplus,alt+residual-boundary+tensor(sb,nl)/4)
        reject('delete sterile tensor source',hplus,alt+residual-boundary-tensor(ls,sr)/4)
        reject('delete covariant boundary term',hplus,alt+residual)
        reject('keep old dipole coefficient after changing placement',combined,alternate+residual-boundary-tensor(ls,rr)/4)
        reject('keep old current coefficient after changing placement',combined,alternate+residual-boundary-current/2)
        reject('ordinary derivative IBP drops curvature divergence',qvp,qvr+qvh+boundaryv)


# Exact Grassmann source expansion, including the charged-lepton weak term.
def gv(a,v):return [sum((v[j]*a[i][j] for j in range(4)),Grass()) for i in range(4)]
def gb(v):return [sum((v[j].dagger()*gamma[0][j][i] for j in range(4)),Grass()) for i in range(4)]
def gc(v):return gv(C,gb(v))
def ga(a,b):return [x+y for x,y in zip(a,b)]
def gs(c,v):return [x*c for x in v]
def gsum(xs):
    z=[Grass() for i in range(4)]
    for v in xs:z=ga(z,v)
    return z
def gd(a,b):return sum((x*y for x,y in zip(a,b)),Grass())
def generators(off):return [Grass({(off+2*s,):1}) for s in range(4)]
na,ns,ne=2,3,2;phi=(Q(2,1),Q(-1,2));hc=(-phi[1],phi[0]);y=rm(na,ns);ye=rm(na,ne);f=rm(na,ns)
l=[[gv(PL,generators(16*al+8*a)) for a in range(2)] for al in range(na)]
bl=[[gb(v) for v in row] for row in l]
nn=[ga(generators(64+8*i),gc(generators(64+8*i))) for i in range(ns)]
nl=[gv(PL,n) for n in nn];nr=[gv(PR,n) for n in nn];bnr=[gb(n) for n in nr]
ee=[gv(PR,generators(128+8*r)) for r in range(ne)];be=[gb(e) for e in ee]
xb=[gsum(gs(phi[a],bl[al][a]) for a in range(2)) for al in range(na)];xc=[gv(C,x) for x in xb]
sr=[gsum(gs(y[al][i],xc[al]) for al in range(na)) for i in range(ns)]
sb=[[gsum([gs(y[al][j].conj()*phi[a].conj(),bnr[j]) for j in range(ns)]+[gs(ye[al][r].conj()*hc[a],be[r]) for r in range(ne)]) for a in range(2)] for al in range(na)]
hn=mul(dagger(y),f);he=mul(dagger(ye),f);kk=mul(f,transpose(y));mass=matrix(ns,ns,lambda i,j:(2,3,5)[i] if i==j else 0);fm=mul(f,mass)
ff=[[zero(2) for n in range(4)] for m in range(4)]
for m,n in pairs:ff[m][n]=herm();ff[n][m]=scale(-1,ff[m][n])
def gt(lrow,z,m,n):return gd(lrow,gv(sigma[m][n],z))
lepton=Grass();sterile=Grass();charged=Grass();nsource=Grass();masspiece=Grass();active=Grass()
for m,n in product(range(4),repeat=2):
    fp=aw(ff[m][n],phi);scalar_n=sum((phi[a].conj()*fp[a] for a in range(2)),Q());scalar_e=sum((hc[a]*fp[a] for a in range(2)),Q())
    for al,i in product(range(na),range(ns)):
        lepton+=sum((gt(sb[al][a],nl[i],m,n)*(f[al][i]*fp[a]) for a in range(2)),Grass())
        nsource+=sum((gt(bl[al][a],sr[i],m,n)*(f[al][i]*fp[a]) for a in range(2)),Grass())
        for j in range(ns):masspiece+=sum((gt(bl[al][a],nr[j],m,n)*(f[al][i]*mass[i][j]*fp[a]) for a in range(2)),Grass())
    sterile+=sum((gt(bnr[j],nl[i],m,n)*(hn[j][i]*scalar_n) for j,i in product(range(ns),repeat=2)),Grass())
    charged+=sum((gt(be[r],nl[i],m,n)*(he[r][i]*scalar_e) for r in range(ne) for i in range(ns)),Grass())
    active+=sum((gt(bl[al][a],xc[beta],m,n)*(kk[al][beta]*fp[a]) for al,beta in product(range(na),repeat=2) for a in range(2)),Grass())
check('lepton tensor source expands into sterile and charged terms',lepton,sterile+charged)
check('sterile Yukawa tensor source retains full active weak tensor',nsource,active)
directmass=Grass()
for m,n in product(range(4),repeat=2):
    fp=aw(ff[m][n],phi)
    directmass+=sum((gt(bl[al][a],nr[i],m,n)*(fm[al][i]*fp[a]) for al in range(na) for i in range(ns) for a in range(2)),Grass())
check('retained mass tensor has ordered FM coefficient',masspiece,directmass)
for m,n in pairs:
    for j,i in product(range(ns),repeat=2):check('sterile Majorana tensor is flavor antisymmetric',gt(bnr[j],nl[i],m,n),-gt(bnr[i],nl[j],m,n))
cs=matrix(na,na,lambda a,b:Q(a+b+1,a*b+2));sel=mul(cs,tuple(tuple(z.conj() for z in row) for row in y));hsel=mul(dagger(y),sel)
check('selected sterile tensor flavor matrix is symmetric',hsel,transpose(hsel))
selectedsterile=Grass()
for m,n in product(range(4),repeat=2):
    fp=aw(ff[m][n],phi);sc=sum((phi[a].conj()*fp[a] for a in range(2)),Q())
    selectedsterile+=sum((gt(bnr[j],nl[i],m,n)*(hsel[j][i]*sc) for j,i in product(range(ns),repeat=2)),Grass())
check('selected neutrino lepton-source tensor cancels by Grassmann symmetry',selectedsterile,Grass())
check('charged hypercharge source vanishes algebraically',sum((hc[a]*phi[a] for a in range(2)),Q()),Q())
reject('charged weak source incorrectly discarded',sum((hc[a]*aw(T[0],phi)[a] for a in range(2)),Q()),Q())
reject('discard charged-lepton tensor in alternative source map',lepton,sterile)


# The optional placement map uses the actual tensor covector and its
# indexed Dirac adjoint, rather than assigning a Majorana half by hand.
xcol=gsum(gs(f[al][0].conj()*phi[b].conj()*ff[m][n][b][a]/4,gv(sigma[m][n],l[al][a])) for al in range(na) for a,b in product(range(2),repeat=2) for m,n in product(range(4),repeat=2))
xrow=gsum(gs(f[al][0]*aw(ff[m][n],phi)[a]/4,[sum((bl[al][a][s]*sigma[m][n][s][t] for s in range(4)),Grass()) for t in range(4)]) for al in range(na) for a in range(2) for m,n in product(range(4),repeat=2))
check('placement map indexed x is the actual tensor row adjoint',gb(xcol),xrow)
check('placement map x is left chiral',gv(PL,xcol),xcol)
z=generators(192);en=ga(z,gc(z));er=gv(PR,en);xreal=ga(xcol,gc(xcol));packet=gd(xrow,er)
check('placement map Majorana completion cancels both channels',packet+packet.dagger(),gd(gb(xreal),en))
for a,b in zip(xreal,gc(xreal)):check('placement map preserves Majorana reality',a,b)
reject('extra half in placement Majorana map',packet+packet.dagger(),gd(gb(xreal),en)*R(1,2))


result={
    'status':'passed','python':sys.version.split()[0],
    'input_candidate':'4ed8560e019dfe74ed6c5178b0f990fe3ccb3008',
    'arithmetic':'Exact Gaussian rational Clifford, Lorentz/weak covariance, connection jets and exterior Grassmann algebra; standard library only.',
    'positive_check_count':sum(counts.values()),'checks_by_category':dict(sorted(counts.items())),
    'negative_control_count':len(mutations),'negative_controls_rejected':sorted(mutations),
    'lorentz_covariance_system':{'unknowns':len(keys),'equations':len(rows),'rank':lorank,'nullity':len(keys)-lorank},
    'weak_covariance_systems':weak_receipts,'circular_Fourier_modes':circular,'raw_F_coefficient_image':source_image_receipts,
    'trial_configuration':{'generic_noncommuting_connection_IBP_source_trials':8,'active_flavors':2,'sterile_flavors':3,'charged_lepton_flavors':2,'complete_triple_gamma_components':64,'complete_Lorentz_and_weak_coefficient_systems':True},
    'scope':'Algebraic derivative-onphi basis and exact sourced placement transport of selected packet; not a complete derivative-EFT quotient or a physical helicity amplitude/Wilson coefficient/beta function.'
}
Path(__file__).with_name('CURVATURE_DERIVATIVE_RESULTS.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps(result,indent=2,sort_keys=True))
