#!/usr/bin/env python3
"""Exact action/source audit for the selected matrix electroweak packet.

Only standard-library imports. No loop formula or repository module is
imported. All failures are explicit exceptions and survive Python -O.
"""
from collections import Counter
from fractions import Fraction as R
from itertools import permutations
import json
from pathlib import Path
import random
import sys


class Q:
    __slots__ = ("a", "b")
    def __init__(self, a=0, b=0):
        if isinstance(a, Q): self.a, self.b = a.a, a.b
        else: self.a, self.b = R(a), R(b)
    def __add__(self, o):
        o = Q(o); return Q(self.a+o.a, self.b+o.b)
    __radd__ = __add__
    def __neg__(self): return Q(-self.a, -self.b)
    def __sub__(self, o): return self + -Q(o)
    def __rsub__(self, o): return Q(o) + -self
    def __mul__(self, o):
        if isinstance(o, Dual): return o*self
        o = Q(o); return Q(self.a*o.a-self.b*o.b, self.a*o.b+self.b*o.a)
    __rmul__ = __mul__
    def __truediv__(self, o):
        o=Q(o); den=o.a*o.a+o.b*o.b
        return Q((self.a*o.a+self.b*o.b)/den,(self.b*o.a-self.a*o.b)/den)
    def __eq__(self, o):
        o=Q(o); return self.a==o.a and self.b==o.b
    def __bool__(self): return bool(self.a or self.b)
    def conj(self): return Q(self.a,-self.b)
    def __repr__(self): return "Q(%s,%s)" % (self.a,self.b)


class Dual:
    """Forward derivative of a literal action in independent real jets."""
    def __init__(self, v=0, index=None, d=None):
        self.v=Q(v); self.d=dict(d or {})
        if index is not None: self.d[index]=Q(1)
    @staticmethod
    def cast(x): return x if isinstance(x,Dual) else Dual(x)
    def __add__(self,o):
        o=Dual.cast(o); d=dict(self.d)
        for k,v in o.d.items(): d[k]=d.get(k,Q())+v
        return Dual(self.v+o.v,d={k:v for k,v in d.items() if v})
    __radd__=__add__
    def __neg__(self): return Dual(-self.v,d={k:-v for k,v in self.d.items()})
    def __sub__(self,o): return self+-Dual.cast(o)
    def __rsub__(self,o): return Dual.cast(o)+-self
    def __mul__(self,o):
        o=Dual.cast(o); d={k:v*o.v for k,v in self.d.items()}
        for k,v in o.d.items(): d[k]=d.get(k,Q())+self.v*v
        return Dual(self.v*o.v,d={k:v for k,v in d.items() if v})
    __rmul__=__mul__
    def conj(self): return Dual(self.v.conj(),d={k:v.conj() for k,v in self.d.items()})


I=Q(0,1); eta=(1,-1,-1,-1)
def mat(n,m,f): return tuple(tuple(Q(f(i,j)) for j in range(m)) for i in range(n))
def zero(n,m=None): return mat(n,n if m is None else m,lambda i,j:0)
def eye(n): return mat(n,n,lambda i,j:i==j)
def add(a,b): return tuple(tuple(x+y for x,y in zip(ar,br)) for ar,br in zip(a,b))
def scale(c,a): return tuple(tuple(c*x for x in ar) for ar in a)
def sub(a,b): return add(a,scale(-1,b))
def mul(a,b):
    return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(len(b)) if a[i][k] and b[k][j]),Q())
                       for j in range(len(b[0]))) for i in range(len(a)))
def dagger(a): return tuple(tuple(a[j][i].conj() for j in range(len(a))) for i in range(len(a[0])))
def comm(a,b): return sub(mul(a,b),mul(b,a))
def msum(xs,n=4,m=None):
    z=zero(n,m)
    for x in xs: z=add(z,x)
    return z
def tr(a): return sum((a[i][i] for i in range(len(a))),Q())
id2,id4=eye(2),eye(4)
pauli=(mat(2,2,lambda i,j: int(i!=j)),((Q(),-I),(I,Q())),mat(2,2,lambda i,j:(1 if i==0 else -1) if i==j else 0))
T=tuple(scale(R(1,2),p) for p in pauli)
g0=mat(4,4,lambda i,j:(1 if i<2 else -1) if i==j else 0)
def spatial(p):
    return tuple(tuple(p[i][j-2] if i<2 and j>=2 else -p[i-2][j] if i>=2 and j<2 else Q() for j in range(4)) for i in range(4))
gamma=(g0,)+tuple(spatial(p) for p in pauli)
sigma=tuple(tuple(scale(I/2,comm(gamma[m],gamma[n])) for n in range(4)) for m in range(4))
g5=scale(I,mul(mul(mul(gamma[0],gamma[1]),gamma[2]),gamma[3]))
PL=scale(R(1,2),sub(id4,g5)); PR=scale(R(1,2),add(id4,g5))
C=scale(I,mul(gamma[2],gamma[0]))
def eps(a,b,c):
    if len({a,b,c})<3: return 0
    return 1 if (a,b,c) in ((0,1,2),(1,2,0),(2,0,1)) else -1
def grade3(m,n,r):
    labels=(m,n,r)
    return msum((scale((-1)**sum(p[a]>p[b] for a in range(3) for b in range(a+1,3)),mul(mul(gamma[labels[p[0]]],gamma[labels[p[1]]]),gamma[labels[p[2]]])) for p in permutations(range(3))),4) if len({m,n,r})==3 else zero(4)
G=tuple(tuple(tuple(scale(R(1,6),grade3(m,n,r)) for r in range(4)) for n in range(4)) for m in range(4))
def slash(v): return msum((scale(v[m],gamma[m]) for m in range(4)))
def adj(a):
    return mul(gamma[0],dagger(a)) if len(a)==1 else mul(dagger(a),gamma[0])
def applyweak(a,phi): return tuple(sum((a[i][j]*phi[j] for j in range(2)),Q()) for i in range(2))
def weakrow(l,a): return tuple(msum((scale(a[b][c],l[b]) for b in range(2)),1,4) for c in range(2))
def bil(l,ph,n,s=id4):
    return sum((mul(mul(l[a],s),n)[0][0]*ph[a] for a in range(2)),Q())
def weakadd(a,b): return tuple(add(a[i],b[i]) for i in range(2))
def weakscale(c,a): return tuple(scale(c,x) for x in a)
def weakmsum(xs):
    z=(zero(1,4),zero(1,4))
    for x in xs:z=weakadd(z,x)
    return z
def vadd(a,b): return tuple(x+y for x,y in zip(a,b))
def vscale(c,a): return tuple(c*x for x in a)
def vsum(xs):
    z=(Q(),Q())
    for x in xs:z=vadd(z,x)
    return z


counts=Counter(); mutations=set(); rng=random.Random(2026101507)
def check(label,a,b=0):
    if a!=b: raise RuntimeError("check failed: "+label+"; "+repr(a)+" != "+repr(b))
    counts[label]+=1
def reject(label,a,b=0):
    if a==b: raise RuntimeError("mutation survived: "+label)
    mutations.add(label)
def real(): return R(rng.randint(-3,3),rng.randint(1,3))
def scalar(): return Q(real(),real())
def rm(n,m): return mat(n,m,lambda i,j:scalar())
def row(P=PR): return mul(rm(1,4),P)
def col(P=PL): return mul(P,rm(4,1))
def herm(): return msum((scale(real(),x) for x in (id2,)+T),2)
def dot4(a,b): return sum((eta[m]*a[m]*b[m] for m in range(4)),Q())


for m in range(4):
    for n in range(4):
        check("Clifford",add(mul(gamma[m],gamma[n]),mul(gamma[n],gamma[m])),scale(2*eta[m] if m==n else 0,id4))
    transpose=tuple(zip(*gamma[m]))
    check("charge conjugation",mul(mul(C,transpose),scale(-1,C)),scale(-1,gamma[m]))
for a in range(3):
    for b in range(3):
        check("weak generator commutator",comm(T[a],T[b]),msum((scale(I*eps(a,b,c),T[c]) for c in range(3)),2))
for a in range(2):
    for b in range(2):
        for c in range(2):
            for d in range(2):
                check("weak generator normalization",sum((t[a][b]*t[c][d] for t in T),Q()),Q(R(1,2)*int(a==d and b==c)-R(1,4)*int(a==b and c==d)))


gp=R(3,5); g=R(2,3); yl=R(-1,2)


for trial in range(6):
    # Differentiate the gauge action in its actual component fields.
    w=[ [Dual(real(),("w",a,m)) for m in range(4)] for a in range(3)]
    dw=[[[Dual(real(),("dw",a,m,n)) for n in range(4)] for m in range(4)] for a in range(3)]
    fw=[[[dw[a][n][m]-dw[a][m][n]-g*sum((eps(b,c,a)*w[b][m]*w[c][n] for b in range(3) for c in range(3)),Dual()) for n in range(4)] for m in range(4)] for a in range(3)]
    lag=-R(1,4)*sum((eta[m]*eta[n]*fw[a][m][n]*fw[a][m][n] for a in range(3) for m in range(4) for n in range(4)),Dual())
    for a in range(3):
        for m in range(4):
            algebra=-g*sum((eps(b,c,a)*w[b][n].v*eta[n]*eta[m]*fw[c][n][m].v for b in range(3) for c in range(3) for n in range(4)),Q())
            check("weak gauge action algebraic Euler term",lag.d.get(("w",a,m),Q()),algebra)
            for n in range(4):
                check("weak gauge action derivative Euler term",lag.d.get(("dw",a,m,n),Q()),-eta[n]*eta[m]*fw[a][n][m].v)
            if trial==0 and a==0 and m==1:
                reject("wrong Dplus adjoint sign",lag.d.get(("w",a,m),Q()),-algebra)
    # Literal scalar kinetic action, with independent real B and W variables.
    b=[Dual(real(),("b",m)) for m in range(4)]
    ph=tuple(scalar() for a in range(2)); dph=tuple(tuple(scalar() for a in range(2)) for m in range(4))
    dc=[]
    for m in range(4):
        am=[[gp*yl*b[m]*(int(a==bb))+g*sum((T[c][a][bb]*w[c][m] for c in range(3)),Dual()) for bb in range(2)] for a in range(2)]
        dc.append([Dual(dph[m][a])+I*sum((am[a][bb]*ph[bb] for bb in range(2)),Dual()) for a in range(2)])
    kinetic=sum((eta[m]*dc[m][a].conj()*dc[m][a] for m in range(4) for a in range(2)),Dual())
    for m in range(4):
        dy=I*yl*eta[m]*sum((ph[a].conj()*dc[m][a].v-dc[m][a].v.conj()*ph[a] for a in range(2)),Q())
        check("scalar hypercharge action current",kinetic.d.get(("b",m),Q()),-gp*dy)
        for a in range(3):
            jw=I*eta[m]*sum((ph[c].conj()*T[a][c][d]*dc[m][d].v-dc[m][c].v.conj()*T[a][c][d]*ph[d] for c in range(2) for d in range(2)),Q())
            check("scalar weak action current",kinetic.d.get(("w",a,m),Q()),-g*jw)
            if trial==0 and m==1 and a==1:
                reject("weak current uses tau instead of tau over two",kinetic.d.get(("w",a,m),Q()),-2*g*jw)
    # The connection part of bar L i slash D L is differentiated before
    # collecting the fermion-current generators; spinor slots stay ordered.
    fl=tuple(row() for a in range(2)); fc=tuple(col() for a in range(2))
    fermion=Dual()
    for m in range(4):
        for a in range(2):
            for bb in range(2):
                connection=gp*yl*b[m]*int(a==bb)+g*sum((T[c][a][bb]*w[c][m] for c in range(3)),Dual())
                fermion-=eta[m]*connection*mul(mul(fl[a],gamma[m]),fc[bb])[0][0]
    for m in range(4):
        fy=yl*eta[m]*sum((mul(mul(fl[a],gamma[m]),fc[a])[0][0] for a in range(2)),Q())
        check("fermion hypercharge action current",fermion.d.get(("b",m),Q()),-gp*fy)
        for a in range(3):
            fwcurrent=eta[m]*sum((T[a][bb][c]*mul(mul(fl[bb],gamma[m]),fc[c])[0][0] for bb in range(2) for c in range(2)),Q())
            check("fermion weak action current",fermion.d.get(("w",a,m),Q()),-g*fwcurrent)


def connection_jets():
    aa=[herm() for m in range(4)]
    da=[[herm() for n in range(4)] for m in range(4)]  # da[field index][derivative]
    dda={}
    for m in range(4):
        for n in range(4):
            for r in range(n,4): dda[m,n,r]=herm()
    dd=lambda m,n,r:dda[m,min(n,r),max(n,r)]
    ff=[[add(sub(da[n][m],da[m][n]),scale(I,comm(aa[m],aa[n]))) for n in range(4)] for m in range(4)]
    df=[[[add(add(sub(dd(n,m,r),dd(m,n,r)),scale(I,add(comm(da[m][r],aa[n]),comm(aa[m],da[n][r])))),scale(I,comm(aa[r],ff[m][n]))) for n in range(4)] for m in range(4)] for r in range(4)]
    return aa,da,ff,df


def matter_jets(aa,da):
    l=tuple(row() for a in range(2)); ph=tuple(scalar() for a in range(2))
    dl=[tuple(row() for a in range(2)) for m in range(4)]
    dp=[tuple(scalar() for a in range(2)) for m in range(4)]
    ll={(m,n):tuple(row() for a in range(2)) for m in range(4) for n in range(m,4)}
    pp={(m,n):tuple(scalar() for a in range(2)) for m in range(4) for n in range(m,4)}
    dlc=[weakadd(dl[m],weakscale(-I,weakrow(l,aa[m]))) for m in range(4)]
    dpc=[vadd(dp[m],vscale(I,applyweak(aa[m],ph))) for m in range(4)]
    dl2=[]; dp2=[]
    for m in range(4):
        lr=[]; pr=[]
        for n in range(4):
            lr.append(weakmsum((ll[min(m,n),max(m,n)],weakscale(-I,weakrow(dl[m],aa[n])),weakscale(-I,weakrow(dl[n],aa[m])),weakscale(-I,weakrow(l,da[n][m])),weakscale(-1,weakrow(l,mul(aa[n],aa[m]))))))
            pr.append(vsum((pp[min(m,n),max(m,n)],vscale(I,applyweak(da[n][m],ph)),vscale(I,applyweak(aa[n],dp[m])),vscale(I,applyweak(aa[m],dp[n])),vscale(-1,applyweak(mul(aa[m],aa[n]),ph)))))
        dl2.append(lr); dp2.append(pr)
    boxl=weakmsum((weakscale(eta[m],dl2[m][m]) for m in range(4)))
    boxp=vsum((vscale(eta[m],dp2[m][m]) for m in range(4)))
    boxlo=weakmsum((weakscale(eta[m],ll[m,m]) for m in range(4)))
    boxpo=vsum((vscale(eta[m],pp[m,m]) for m in range(4)))
    kl=tuple(msum((scale(-I,mul(dlc[m][a],gamma[m])) for m in range(4)),1,4) for a in range(2))
    k2l=tuple(msum((scale(-1,mul(mul(dl2[m][n][a],gamma[n]),gamma[m])) for m in range(4) for n in range(4)),1,4) for a in range(2))
    return dict(l=l,ph=ph,dl=dl,dp=dp,dlc=dlc,dpc=dpc,boxl=boxl,boxp=boxp,boxlo=boxlo,boxpo=boxpo,kl=kl,k2l=k2l)


for trial in range(8):
    aa,da,ff,df=connection_jets()
    for r in range(4):
        for m in range(4):
            for n in range(4):
                check("covariant nonabelian Bianchi jet",add(add(df[r][m][n],df[m][n][r]),df[n][r][m]),zero(2))
    if trial==0:
        missing=msum((scale(I,comm(aa[r],ff[m][n])) for r,m,n in ((0,1,2),(1,2,0),(2,0,1))),2)
        reject("ordinary instead of covariant Bianchi",missing,zero(2))
        reject("commuting weak connection",comm(aa[0],aa[1]),zero(2))
    jet=matter_jets(aa,da)
    l,ph,dlc,dpc=jet['l'],jet['ph'],jet['dlc'],jet['dpc']
    for a in range(2):
        curvature=msum((scale(R(1,2)*ff[m][n][b][a],mul(l[b],sigma[m][n])) for m in range(4) for n in range(4) for b in range(2)),1,4)
        check("matrix row Dirac square",jet['k2l'][a],sub(scale(-1,jet['boxl'][a]),curvature))
        if trial==0 and a==0: reject("opposite row curvature sign",jet['k2l'][a],add(scale(-1,jet['boxl'][a]),curvature))
    n=col(); dn=[rm(4,1) for m in range(4)]
    an=msum((scale(I,mul(gamma[m],dn[m])) for m in range(4)),4,1)
    qsig=Q(); qdiv=Q(); qv=Q(); qe=Q(); qan=Q(); totaldiv=Q()
    for m in range(4):
        for nidx in range(4):
            fp=applyweak(ff[m][nidx],ph)
            qan+=bil(l,fp,an,sigma[m][nidx])
            qdiv+=eta[nidx]*bil(l,applyweak(df[nidx][nidx][m],ph),n,gamma[m])
            qv-=eta[nidx]*bil(l,applyweak(ff[m][nidx],dpc[nidx]),n,gamma[m])
            for r in range(4):
                sg=mul(sigma[m][nidx],gamma[r])
                first=bil(dlc[r],fp,n,sg)
                qsig+=I*R(1,4)*first
                qe-=R(1,2)*bil(l,applyweak(ff[m][nidx],dpc[r]),n,G[m][nidx][r])
                totaldiv+=I*R(1,4)*(first+bil(l,applyweak(df[r][m][nidx],ph),n,sg)+bil(l,applyweak(ff[m][nidx],dpc[r]),n,sg)+bil(l,fp,dn[r],sg))
    od=qsig+qdiv/3
    red=-(qv+qe)/2-qdiv/6-qan/4
    check("ordered matrix operator IBP from product derivative",od,red+totaldiv)
    if trial==0:
        reject("isolated tensor coefficient sign reversed",od,red+totaldiv+qan/2)
        reject("matrix divergence sign reversed",od,red+totaldiv+qdiv/3)
        # Reverse only the actual insertion order F Dphi -> Dphi with F^T.
        wrong=Q()
        for m in range(4):
            for nidx in range(4):
                ftrans=tuple(tuple(ff[m][nidx][b][a] for b in range(2)) for a in range(2))
                wrong-=eta[nidx]*bil(l,applyweak(ftrans,dpc[nidx]),n,gamma[m])
        reject("transpose weak curvature without conjugation",qv,wrong)


for trial in range(6):
    # Entire sourced residual, with independent nonzero sources and all
    # sterile iterates; no free field equation is imposed.
    aa,da,ff,df=connection_jets()
    jets=[matter_jets(aa,da) for a in range(2)]
    ph=jets[0]['ph']
    # All active rows share one scalar field and its derivatives.
    for jet in jets[1:]:
        for key in ('ph','dp','dpc','boxp','boxpo'):jet[key]=jets[0][key]
    f=rm(2,3); mass=((Q(2),Q(1),Q(-1)),(Q(1),Q(4),Q(2)),(Q(-1),Q(2),Q(6))) if trial%2 else mat(3,3,lambda i,j:(2,3,5)[i] if i==j else 0)
    mm=mul(mass,mass)
    def cols():return [rm(4,1) for i in range(3)]
    nl,nr,anl,anr,a2nl,a2nr,a3nl,sr,sn,asr,asn,a2sr=[cols() for k in range(12)]
    def massact(z):return [msum((scale(mass[i][j],z[j]) for j in range(3)),4,1) for i in range(3)]
    def plus(a,b):return [add(x,y) for x,y in zip(a,b)]
    def minus(a,b):return [sub(x,y) for x,y in zip(a,b)]
    er=minus(minus(anl,massact(nr)),sr); el=minus(minus(anr,massact(nl)),sn)
    aer=minus(minus(a2nl,massact(anr)),asr); ael=minus(minus(a2nr,massact(anl)),asn)
    a2er=minus(minus(a3nl,massact(a2nr)),a2sr)
    rr=plus(massact(nr),sr)
    tt=plus(plus(massact(massact(nl)),massact(sn)),asr)
    uu=plus(plus(plus(massact(massact(massact(nr))),massact(massact(sr))),massact(asn)),a2sr)
    for i in range(3):
        check("sterile first iterate residual",sub(anl[i],rr[i]),er[i])
        check("sterile second iterate residual",sub(a2nl[i],tt[i]),add(massact(el)[i],aer[i]))
        check("sterile third iterate residual",sub(a3nl[i],uu[i]),add(add(massact(massact(er))[i],massact(ael)[i]),a2er[i]))
    sbar=[tuple(row(PL) for a in range(2)) for al in range(2)]
    ksbar=[tuple(row() for a in range(2)) for al in range(2)]
    sph=tuple(scalar() for a in range(2)); m2=real(); lam=real(); rho=sum((p.conj()*p for p in ph),Q()); cv=2*m2/3+lam*rho/3
    ebar=[tuple(sub(jets[al]['kl'][a],sbar[al][a]) for a in range(2)) for al in range(2)]
    kebar=[tuple(sub(jets[al]['k2l'][a],ksbar[al][a]) for a in range(2)) for al in range(2)]
    eph=vadd(vadd(vscale(-1,jets[0]['boxp']),vscale(-(m2+2*lam*rho),ph)),sph)
    def b(rows,vec,z,spin=id4):
        return sum((f[al][i]*bil(rows[al],vec,z[i],spin) for al in range(2) for i in range(3)),Q())
    ll=[j['l'] for j in jets]; kl=[j['kl'] for j in jets]; lapl=[weakscale(-1,j['boxl']) for j in jets]
    def tensor(z):
        return sum((b(ll,applyweak(ff[m][n],ph),z,sigma[m][n]) for m in range(4) for n in range(4)),Q())
    jy=[real() for m in range(4)]; jw=[[real() for m in range(4)] for a in range(3)]
    div=[msum((scale(eta[n]*eta[m],df[n][n][m]) for n in range(4)),2) for m in range(4)]
    divb=[tr(v)/(2*gp*yl) for v in div]
    divw=[[2*tr(mul(T[a],div[m]))/g for m in range(4)] for a in range(3)]
    eb=[divb[m]-gp*jy[m] for m in range(4)]; ew=[[divw[a][m]-g*jw[a][m] for m in range(4)] for a in range(3)]
    current=Q(); weak_current=Q(); gauge_res=Q(); literal_div=Q()
    for m in range(4):
        j=eta[m]*b(ll,ph,nl,gamma[m]); ji=[eta[m]*b(ll,applyweak(T[a],ph),nl,gamma[m]) for a in range(3)]
        current-=gp*gp*yl*j*jy[m]/6+g*g*sum((ji[a]*jw[a][m] for a in range(3)),Q())/6
        weak_current-=g*g*sum((ji[a]*jw[a][m] for a in range(3)),Q())/6
        gauge_res-=gp*yl*j*eb[m]/6+g*sum((ji[a]*ew[a][m] for a in range(3)),Q())/6
        literal_div-=eta[m]*b(ll,applyweak(div[m],ph),nl,gamma[m])/6
        split=add(scale(gp*yl*(eb[m]+gp*jy[m]),id2),msum((scale(g*(ew[a][m]+g*jw[a][m]),T[a]) for a in range(3)),2))
        check("matrix gauge divergence source split",div[m],split)
    # Common Higgs-derivative packet cancels identically, so omit it from
    # both sides here; its raw action IBP was checked above independently.
    original=R(5,6)*b(ll,ph,a3nl)+b(lapl,ph,anl)/3+b(ll,vscale(-1,jets[0]['boxp']),anl)/6-b(kl,ph,a2nl)/6+m2*b(ll,ph,anl)/2+m2*b(kl,ph,nl)/2+literal_div-tensor(anl)/4
    reduced=R(5,6)*b(ll,ph,uu)+cv*b(ll,ph,rr)+b(ksbar,ph,rr)/3-b(ll,sph,rr)/6-b(sbar,ph,tt)/6+m2*b(sbar,ph,nl)/2+current-tensor(rr)/12
    residual=b(kebar,ph,anl)/3+b(ll,eph,anl)/6-b(ebar,ph,a2nl)/6+m2*b(ebar,ph,nl)/2+gauge_res+R(5,6)*b(ll,ph,plus(plus(massact(massact(er)),massact(ael)),a2er))+cv*b(ll,ph,er)+b(ksbar,ph,er)/3-b(ll,sph,er)/6-b(sbar,ph,plus(massact(el),aer))/6-tensor(er)/12
    check("complete matrix sourced residual certificate",original,reduced+residual)
    if trial==0:
        reject("combined dipole replaced by isolated minus quarter",original,reduced+residual-tensor(rr)/6)
        reject("sterile tensor source omitted",original,reduced+residual+tensor(sr)/12)
        reject("weak and hypercharge current packet omitted",original,reduced+residual-current)
        reject("weak current coefficient halved",original,reduced+residual-weak_current/2)
        reject("scalar source omitted",original,reduced+residual+b(ll,sph,rr)/6)
        reject("differentiated lepton source omitted",original,reduced+residual-b(ksbar,ph,rr)/3)


for trial in range(4):
    aa,da,ff,df=connection_jets(); jets=[matter_jets(aa,da) for al in range(2)]
    ph=jets[0]['ph']; dp=jets[0]['dp']; dpc=jets[0]['dpc']; bp=jets[0]['boxp']; bpo=jets[0]['boxpo']
    f=rm(2,3); mass=((Q(2),Q(1),Q(-1)),(Q(1),Q(4),Q(2)),(Q(-1),Q(2),Q(6)))
    fm=mul(f,mass); fm2=mul(fm,mass); cv=real(); sph=tuple(scalar() for a in range(2))
    sb=[tuple(row(PL) for a in range(2)) for al in range(2)]
    ds=[[tuple(row(PL) for a in range(2)) for m in range(4)] for al in range(2)]
    dsc=[[weakadd(ds[al][m],weakscale(-I,weakrow(sb[al],aa[m]))) for m in range(4)] for al in range(2)]
    for i in range(3):
        w=msum((scale(f[al][i]*ph[a],jets[al]['l'][a]) for al in range(2) for a in range(2)),1,4)
        ww=msum((scale(fm2[al][i]*ph[a],jets[al]['l'][a]) for al in range(2) for a in range(2)),1,4)
        v=msum((scale(f[al][i]*ph[a],sb[al][a]) for al in range(2) for a in range(2)),1,4)
        kw=msum((scale(-I*fm[al][i],mul(add(scale(ph[a],jets[al]['dl'][m][a]),scale(dp[m][a],jets[al]['l'][a])),gamma[m])) for al in range(2) for a in range(2) for m in range(4)),1,4)
        k2w=msum((scale(-f[al][i],add(add(scale(ph[a],jets[al]['boxlo'][a]),msum((scale(2*eta[m]*dp[m][a],jets[al]['dl'][m][a]) for m in range(4)),1,4)),scale(bpo[a],jets[al]['l'][a]))) for al in range(2) for a in range(2)),1,4)
        kv=msum((scale(-I*f[al][i],mul(add(scale(ph[a],ds[al][m][a]),scale(dp[m][a],sb[al][a])),gamma[m])) for al in range(2) for a in range(2) for m in range(4)),1,4)
        ks=msum((scale(-I*f[al][i]*ph[a],mul(dsc[al][m][a],gamma[m])) for al in range(2) for a in range(2) for m in range(4)),1,4)
        sp=msum((scale(f[al][i]*sph[a],jets[al]['l'][a]) for al in range(2) for a in range(2)),1,4)
        tf=msum((scale(f[al][i]*ff[m][n][a][b]*ph[b],mul(jets[al]['l'][a],sigma[m][n])) for al in range(2) for a in range(2) for b in range(2) for m in range(4) for n in range(4)),1,4)
        xr=msum((scale(R(5,6),add(ww,k2w)),scale(cv,w),scale(R(1,3),ks),scale(R(-1,6),sp),scale(R(-1,6),kv),scale(R(-1,12),tf)),1,4)
        yr=sub(scale(R(5,6),kw),msum((scale(fm[al][i]*ph[a]/6,sb[al][a]) for al in range(2) for a in range(2)),1,4))
        # Independent explicitly indexed columns (8), using full matrix D.
        xc=msum((scale(R(5,6)*fm2[al][i].conj()*ph[a].conj()+cv*f[al][i].conj()*ph[a].conj(),adj(jets[al]['l'][a])) for al in range(2) for a in range(2)),4,1)
        xc=add(xc,msum((scale(-R(5,6)*f[al][i].conj(),add(add(scale(ph[a].conj(),adj(jets[al]['boxl'][a])),msum((scale(2*eta[m]*dpc[m][a].conj(),adj(jets[al]['dlc'][m][a])) for m in range(4)),4,1)),scale(bp[a].conj(),adj(jets[al]['l'][a])))) for al in range(2) for a in range(2)),4,1))
        xc=add(xc,msum((scale(I*f[al][i].conj()/6,mul(gamma[m],sub(scale(ph[a].conj(),adj(dsc[al][m][a])),scale(dpc[m][a].conj(),adj(sb[al][a]))))) for al in range(2) for a in range(2) for m in range(4)),4,1))
        xc=sub(xc,msum((scale(f[al][i].conj()*sph[a].conj()/6,adj(jets[al]['l'][a])) for al in range(2) for a in range(2)),4,1))
        tcol=msum((scale(f[al][i].conj()*ph[b].conj()*ff[m][n][b][a]/12,mul(sigma[m][n],adj(jets[al]['l'][a]))) for al in range(2) for a in range(2) for b in range(2) for m in range(4) for n in range(4)),4,1)
        xc=sub(xc,tcol)
        yc=sub(msum((scale(5*I*fm[al][i].conj()/6,mul(gamma[m],add(scale(ph[a].conj(),adj(jets[al]['dlc'][m][a])),scale(dpc[m][a].conj(),adj(jets[al]['l'][a]))))) for al in range(2) for a in range(2) for m in range(4)),4,1),msum((scale(fm[al][i].conj()*ph[a].conj()/6,adj(sb[al][a])) for al in range(2) for a in range(2)),4,1))
        check("fully indexed matrix x column adjoint",xc,adj(xr))
        check("fully indexed matrix y column adjoint",yc,adj(yr))
        check("x left chirality",mul(PL,xc),xc); check("y right chirality",mul(PR,yc),yc)
        if trial==0 and i==0:
            wrong=msum((scale(f[al][i].conj()*ph[b].conj()*ff[m][n][a][b]/12,mul(sigma[m][n],adj(jets[al]['l'][a]))) for al in range(2) for a in range(2) for b in range(2) for m in range(4) for n in range(4)),4,1)
            reject("x curvature adjoint uses ab instead of ba",xc,add(xc,sub(tcol,wrong)))
            reject("Majorana sterile source derivative lost",xr,add(xr,scale(R(1,6),kv)))
            wrong_y=msum((scale(-5*I*(f[al][i]*sum((mass[j][i] for j in range(3)),Q()))/6,mul(add(scale(ph[a],jets[al]['dl'][m][a]),scale(dp[m][a],jets[al]['l'][a])),gamma[m])) for al in range(2) for a in range(2) for m in range(4)),1,4)
            wrong_y=sub(wrong_y,msum((scale(f[al][i]*sum((mass[j][i] for j in range(3)),Q())*ph[a]/6,sb[al][a]) for al in range(2) for a in range(2)),1,4))
            reject("mass summed after fixing sterile flavor slot",yr,wrong_y)
    # Light-field IBP with the full noncommuting matrix connection.
    h=[real() for m in range(4)]; p=[real() for m in range(4)]; q=[-h[m]-p[m] for m in range(4)]
    erow=tuple(row() for a in range(2)); nl=col(); anl=scale(-1,mul(slash(p),nl)); a2=mul(mul(slash(p),slash(p)),nl); m2=real()
    ker=tuple(msum((mul(sub(scale(q[m],erow[a]),weakrow(erow,aa[m])[a]),gamma[m]) for m in range(4)),1,4) for a in range(2))
    original=bil(ker,ph,anl)/3-bil(erow,ph,a2)/6+m2*bil(erow,ph,nl)/2
    dpf=[vscale(I,vadd(vscale(h[m],ph),applyweak(aa[m],ph))) for m in range(4)]
    ell=[msum((msum((scale(I*dpf[m][a]/3,mul(gamma[m],anl)) for m in range(4)),4,1),scale(ph[a]/6,a2),scale(m2*ph[a]/2,nl)),4,1) for a in range(2)]
    transported=sum((mul(erow[a],ell[a])[0][0] for a in range(2)),Q())
    check("full weak lepton integrated adjoint",original,transported)


for trial in range(6):
    # Conjugation rotates a weak adjoint contact; the connection difference
    # must transform homogeneously even though the connection does not.
    U=scale(R(1,5),sub(scale(3,id2),scale(4*I,pauli[trial%3])))
    check("weak finite transformation unitary",mul(U,dagger(U)),id2)
    ph=tuple(scalar() for a in range(2)); l=tuple(row() for a in range(2)); n=col()
    lp=weakrow(l,dagger(U)); pp=applyweak(U,ph)
    j=[bil(l,applyweak(t,ph),n,gamma[trial%4]) for t in T]
    jp=[bil(lp,applyweak(t,pp),n,gamma[trial%4]) for t in T]
    rot=[[2*tr(mul(mul(mul(T[a],U),T[b]),dagger(U))) for b in range(3)] for a in range(3)]
    for a in range(3):
        check("weak contact adjoint covariance",jp[a],sum((rot[a][b]*j[b] for b in range(3)),Q()))
        shift=g*(j[a]+j[a].conj())/6
        check("weak field shift reality",shift,shift.conj())
        check("weak real shift cancels residual and conjugate",-(g*j[a]/6)-(g*j[a].conj()/6)+shift,Q())
        if trial==0 and a==0:
            reject("half-sized real weak field shift",-(g*j[a]/6)-(g*j[a].conj()/6)+shift/2,Q())
        for b in range(3):
            check("weak contact adjoint real rotation",rot[a][b],rot[a][b].conj())
    invariant=bil(lp,pp,n); check("hypercharge contact weak invariance",invariant,bil(l,ph,n))


class Grass:
    """Finite exterior algebra, with dagger reversing the generator order."""
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
                sign=(-1)**sum(i>j for i in a for j in b);key=tuple(sorted(a+b))
                d[key]=d.get(key,Q())+sign*x*y
        return Grass(d)
    __rmul__=__mul__
    def dagger(self):
        z=Grass()
        for k,v in self.d.items():
            term=Grass.cast(v.conj())
            for j in reversed(k):term=term*Grass({(j^1,):1})
            z+=term
        return z
    def __eq__(self,o):return self.d==Grass.cast(o).d
    def __repr__(self):return repr(self.d)
def gv(a,v):return [sum((v[j]*a[i][j] for j in range(4)),Grass()) for i in range(4)]
def gb(v):return [sum((v[j].dagger()*gamma[0][j][i] for j in range(4)),Grass()) for i in range(4)]
def gc(v):return gv(C,gb(v))
def gadd(a,b):return [x+y for x,y in zip(a,b)]
def gdot(a,b):return sum((x*y for x,y in zip(a,b)),Grass())
w=[Grass({(2*i,):1}) for i in range(4)]; z=[Grass({(8+2*i,):1}) for i in range(4)]
e=gadd(z,gc(z)); wr=gadd(w,gc(w))
for a,b in zip(wr,gc(wr)):check("Majorana completed shift reality",a,b)
packet=gdot(gb(w),e); both=packet+packet.dagger(); actual=gdot(gb(wr),e)
check("Grassmann Majorana residual completion",both,actual)
reject("extra Majorana half factor",both,actual*R(1,2))
check("Majorana two varied mass slots",(gdot(gb(wr),e)+gdot(gb(e),wr))*R(1,2),actual)
for m in range(4):
    first=gdot(gb(wr),gv(scale(I,gamma[m]),e))
    second=-gdot(gb(e),gv(scale(I,gamma[m]),wr))
    check("Majorana two varied kinetic slots after IBP",(first+second)*R(1,2),first)


for m in range(4):
    for n in range(m+1,4):
        # Transposition here does not conjugate: this is the coefficient of
        # two odd spinor slots in sigma C, not a Hermitian bilinear.
        spin=mul(mul(mul(PR,sigma[m][n]),C),tuple(zip(*PR)))
        check("tensor Grassmann combined-slot symmetry",tuple(zip(*spin)),spin)
ph=(Q(3,0),Q())
for t in (id2,)+T:
    tp=applyweak(t,ph)
    full=tuple(tuple(tp[a]*ph[b] for b in range(2)) for a in range(2))
    anti=tuple(tuple((full[a][b]-full[b][a])/2 for b in range(2)) for a in range(2))
    check("weak source combined-slot antisymmetry",tuple(zip(*anti)),scale(-1,anti))
    if t==id2:check("single-flavor hypercharge source vanishes",anti,zero(2))
    if t==T[0]:reject("flavor-only antisymmetrization erases weak source",anti,zero(2))


result={
    "status":"passed",
    "source_candidate":"d822e9f4d70ad5b455436f40abc3c2b45909bbef",
    "python":sys.version.split()[0],
    "arithmetic":"exact Gaussian rational arithmetic; standard library only",
    "checks_by_category":dict(sorted(counts.items())),
    "positive_check_count":sum(counts.values()),
    "negative_controls_rejected":sorted(mutations),
    "negative_control_count":len(mutations),
    "trial_configuration":{"gauge_and_scalar_action_trials":6,"noncommuting_connection_and_IBP_trials":8,"full_sourced_residual_trials":6,"indexed_adjoint_trials":4,"active_flavors":2,"sterile_flavors":3,"weak_components":2,"finite_weak_covariance_trials":6},
    "scope":"Action algebra for the complete selected matched packet; no new loop integration, full EFT coefficient map, finite amplitude or physical running.",
}
Path(__file__).with_name("ELECTROWEAK_REDUCTION_RESULTS.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
print(json.dumps(result,indent=2,sort_keys=True))
