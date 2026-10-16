#!/usr/bin/env python3
"""Exact scalar-current tensor, IBP, sourced-EOM and flavor audit.

Standard library only; complete weak polynomial covariance systems,
Gaussian rational connection jets and finite exterior Grassmann algebra.
No repository imports, loop integration or sampled-only basis argument.
"""
from collections import Counter
from fractions import Fraction as R
from itertools import product
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
        if isinstance(o,(Poly,Grass)):return o*self
        o=Q(o);return Q(self.a*o.a-self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=Q(o);d=o.a*o.a+o.b*o.b
        return Q((self.a*o.a+self.b*o.b)/d,(self.b*o.a-self.a*o.b)/d)
    def __eq__(self,o):o=Q(o);return self.a==o.a and self.b==o.b
    def __bool__(self):return bool(self.a or self.b)
    def conj(self):return Q(self.a,-self.b)
    def __repr__(self):return 'Q(%s,%s)'%(self.a,self.b)


class Poly:
    def __init__(self,d=None):self.d={k:Q(v) for k,v in (d or {}).items() if v}
    @staticmethod
    def cast(x):return x if isinstance(x,Poly) else Poly({():Q(x)})
    @staticmethod
    def var(i):return Poly({(i,):1})
    def __add__(self,o):
        o=Poly.cast(o);d=dict(self.d)
        for k,v in o.d.items():d[k]=d.get(k,Q())+v
        return Poly(d)
    __radd__=__add__
    def __neg__(self):return Poly({k:-v for k,v in self.d.items()})
    def __sub__(self,o):return self+-Poly.cast(o)
    def __mul__(self,o):
        o=Poly.cast(o);d={}
        for a,x in self.d.items():
            for b,y in o.d.items():
                k=tuple(sorted(a+b));d[k]=d.get(k,Q())+x*y
        return Poly(d)
    __rmul__=__mul__
    def __eq__(self,o):return self.d==Poly.cast(o).d
    def __repr__(self):return repr(self.d)


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
                sign=(-1)**sum(i>j for i in a for j in b);k=tuple(sorted(a+b))
                d[k]=d.get(k,Q())+sign*x*y
        return Grass(d)
    __rmul__=__mul__
    def dagger(self):
        out=Grass()
        for k,v in self.d.items():
            t=Grass.cast(v.conj())
            for j in reversed(k):t=t*Grass({(j^1,):1})
            out+=t
        return out
    def __eq__(self,o):return self.d==Grass.cast(o).d
    def __repr__(self):return repr(self.d)


I=Q(0,1);eta=(1,-1,-1,-1)
def mat(n,m,f):return tuple(tuple(Q(f(a,b)) for b in range(m)) for a in range(n))
def eye(n):return mat(n,n,lambda a,b:a==b)
def zero(n,m=None):return mat(n,n if m is None else m,lambda a,b:0)
def add(a,b):return tuple(tuple(x+y for x,y in zip(ar,br)) for ar,br in zip(a,b))
def scale(c,a):return tuple(tuple(c*x for x in ar) for ar in a)
def sub(a,b):return add(a,scale(-1,b))
def mul(a,b):return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(len(b))),Q()) for j in range(len(b[0]))) for i in range(len(a)))
def dagger(a):return tuple(tuple(a[j][i].conj() for j in range(len(a))) for i in range(len(a[0])))
def msum(xs,n=4,m=None):
    z=zero(n,m)
    for x in xs:z=add(z,x)
    return z
pauli=(mat(2,2,lambda a,b:a!=b),((Q(),-I),(I,Q())),((Q(1),Q()),(Q(),Q(-1))))
T=tuple(scale(R(1,2),p) for p in pauli)
g0=mat(4,4,lambda a,b:(1 if a<2 else -1) if a==b else 0)
def spatial(p):return tuple(tuple(p[a][b-2] if a<2 and b>=2 else -p[a-2][b] if a>=2 and b<2 else Q() for b in range(4)) for a in range(4))
gamma=(g0,)+tuple(spatial(p) for p in pauli)
g5=scale(I,mul(mul(mul(gamma[0],gamma[1]),gamma[2]),gamma[3]))
PL=scale(R(1,2),sub(eye(4),g5));PR=scale(R(1,2),add(eye(4),g5));C=scale(I,mul(gamma[2],gamma[0]))
def adj(a):return mul(gamma[0],dagger(a)) if len(a)==1 else mul(dagger(a),gamma[0])
def apply(a,v):return tuple(sum((a[i][j]*v[j] for j in range(2)),Q()) for i in range(2))
def rowweak(l,a):return tuple(msum((scale(a[b][c],l[b]) for b in range(2)),1,4) for c in range(2))
def bil(l,ph,n,spin=eye(4)):return sum((ph[a]*mul(mul(l[a],spin),n)[0][0] for a in range(2)),Q())


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
        pick=next((i for i in range(r,len(a)) if a[i][c]),None)
        if pick is None:continue
        a[r],a[pick]=a[pick],a[r];v=a[r][c];a[r]=[x/v for x in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][c]:
                f=a[i][c];a[i]=[x-f*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r


# Weak tensor completeness: all polynomial coefficients, no field sampling.
ph=tuple(Poly.var(a) for a in range(2));pc=tuple(Poly.var(2+a) for a in range(2))
dp=tuple(Poly.var(4+a) for a in range(2));dpc=tuple(Poly.var(6+a) for a in range(2))
rho=sum((pc[a]*ph[a] for a in range(2)),Poly())
xp=sum((pc[a]*dp[a] for a in range(2)),Poly());xm=sum((dpc[a]*ph[a] for a in range(2)),Poly())
vA=tuple(rho*dp[a] for a in range(2));vB=tuple(ph[a]*xp for a in range(2));vD=tuple(ph[a]*xm for a in range(2))
for a in range(2):
    raw=Poly()
    for t in T:
        tp=sum((t[a][b]*ph[b] for b in range(2)),Poly())
        current=I*sum((pc[b]*t[b][c]*dp[c]-dpc[b]*t[b][c]*ph[c] for b,c in product(range(2),repeat=2)),Poly())
        raw+=tp*current
    check('complete weak scalar current polynomial',raw,I*R(1,2)*vA[a]-I*R(1,4)*(vB[a]+vD[a]))
    check('complete scalar weak coefficient polynomial',raw*R(-1,6),-I*R(1,12)*vA[a]+I*R(1,24)*(vB[a]+vD[a]))
    hyper=I*R(-1,24)*ph[a]*(xp-xm)
    check('complete scalar hypercharge coefficient polynomial',hyper,-I*R(1,24)*vB[a]+I*R(1,24)*vD[a])
    if a==0:
        reject('discard norm derivative in weak scalar current',raw,I*R(1,2)*vA[a])
        reject('hypercharge scalar current uses sum rather than difference',hyper,-I*R(1,24)*ph[a]*(xp+xm))
check('charged-lepton orthogonality polynomial',-ph[1]*ph[0]+ph[0]*ph[1],Poly())


def covsys(conjugate_derivative=False):
    mons=tuple(sorted(tuple(sorted((a,b,c))) for a,b,c in (product(range(2),range(2,4),range(4,6)) if not conjugate_derivative else ((a,b,c) for a in range(2) for b in range(a,2) for c in range(6,8)))))
    keys=tuple((a,m) for a in range(2) for m in mons);idx={k:i for i,k in enumerate(keys)}
    rows=[[Q() for k in keys] for x in range(3*len(keys))]
    for col,(out,m) in enumerate(keys):
        for gen,t in enumerate(T):
            for slot,old in enumerate(m):
                block=old//2;comp=old%2
                for new in range(2):
                    coeff=t[comp][new] if block in (0,2) else -t[new][comp]
                    mm=list(m);mm[slot]=2*block+new;mm=tuple(sorted(mm))
                    rows[gen*len(keys)+idx[out,mm]][col]+=coeff
            for a in range(2):rows[gen*len(keys)+idx[a,m]][col]-=t[a][out]
    return keys,rows
weak_receipts=[]
for conj,expected,vs in ((False,2,(vA,vB)),(True,1,(vD,))):
    keys,rows=covsys(conj);rr=rank(rows,len(keys))
    check('complete scalar derivative weak covariance nullity',len(keys)-rr,expected)
    cand=[]
    for v in vs:
        coeff=tuple(v[a].d.get(m,Q()) for a,m in keys);cand.append(coeff)
        check('explicit scalar derivative invariant solves covariance',tuple(sum((x*y for x,y in zip(row,coeff)),Q()) for row in rows),tuple(Q() for row in rows))
    check('scalar derivative invariant span is complete',rank(cand,len(keys)),expected)
    weak_receipts.append(dict(derivative='phi_dagger' if conj else 'phi',unknowns=len(keys),equations=len(rows),rank=rr,nullity=len(keys)-rr))
combinedkeys=sorted(set((a,m) for v in (vA,vB,vD) for a in range(2) for m in v[a].d))
check('three scalar derivative contractions are algebraically independent',rank([[v[a].d.get(m,Q()) for a,m in combinedkeys] for v in (vA,vB,vD)],len(combinedkeys)),3)


rng=random.Random(2026101603)
def real():return R(rng.randint(-3,3),rng.randint(1,3))
def scalar():return Q(real(),real())
def rm(n,m):return mat(n,m,lambda a,b:scalar())
def row(P=PR):return mul(rm(1,4),P)
def col(P=PL):return mul(P,rm(4,1))
def herm():return msum((scale(real(),p) for p in (eye(2),)+T),2)


for trial in range(8):
    na,ns=2,3;f=rm(na,ns);m=mat(ns,ns,lambda i,j:(2,3,5)[i] if i==j else 0)
    phi=tuple(scalar() for a in range(2));dphi=[tuple(scalar() for a in range(2)) for mu in range(4)]
    aa=[herm() for mu in range(4)]
    ls=[tuple(row() for a in range(2)) for al in range(na)]
    dl=[[tuple(row() for a in range(2)) for mu in range(4)] for al in range(na)]
    nl=[col() for i in range(ns)];nr=[col(PR) for i in range(ns)]
    dn=[[col() for mu in range(4)] for i in range(ns)]
    sb=[tuple(row(PL) for a in range(2)) for al in range(na)];sr=[col(PR) for i in range(ns)]
    dpcov=[tuple(dphi[mu][a]+I*apply(aa[mu],phi)[a] for a in range(2)) for mu in range(4)]
    dlcov=[[tuple(sub(dl[al][mu][a],scale(I,rowweak(ls[al],aa[mu])[a])) for a in range(2)) for mu in range(4)] for al in range(na)]
    kl=[tuple(msum((scale(-I,mul(dlcov[al][mu][a],gamma[mu])) for mu in range(4)),1,4) for a in range(2)) for al in range(na)]
    an=[msum((scale(I,mul(gamma[mu],dn[i][mu])) for mu in range(4)),4,1) for i in range(ns)]
    er=[sub(sub(an[i],msum((scale(m[i][j],nr[j]) for j in range(ns)),4,1)),sr[i]) for i in range(ns)]
    eb=[tuple(sub(kl[al][a],sb[al][a]) for a in range(2)) for al in range(na)]
    rr=sum((p.conj()*p for p in phi),Q());dr=[sum((dphi[mu][a].conj()*phi[a]+phi[a].conj()*dphi[mu][a] for a in range(2)),Q()) for mu in range(4)]
    def b(rows,phh,nn,sp=eye(4)):return sum((f[al][i]*bil(rows[al],phh,nn[i],sp) for al in range(na) for i in range(ns)),Q())
    j=[b(ls,phi,nl,gamma[mu]) for mu in range(4)]
    direct_div=Q();cov_div=Q();scalar_div=Q();aop=Q();bop=Q();dop=Q();literal=Q()
    g2=R(4,9);gp2=R(9,25);d=(g2+gp2)/24
    for mu in range(4):
        first=sum((f[al][i]*bil(dl[al][mu],phi,nl[i],gamma[mu]) for al in range(na) for i in range(ns)),Q())
        second=b(ls,dphi[mu],nl,gamma[mu])
        third=sum((f[al][i]*bil(ls[al],phi,dn[i][mu],gamma[mu]) for al in range(na) for i in range(ns)),Q())
        direct_div+=first+second+third
        scalar_div+=b(ls,dpcov[mu],nl,gamma[mu])
        aop+=rr*b(ls,dpcov[mu],nl,gamma[mu])
        xpval=sum((phi[a].conj()*dpcov[mu][a] for a in range(2)),Q())
        xmval=sum((dpcov[mu][a].conj()*phi[a] for a in range(2)),Q())
        bop+=j[mu]*xpval;dop+=j[mu]*xmval
        literal-=gp2*R(1,4)*I*j[mu]*(xpval-xmval)/6
        for t in T:
            ji=b(ls,apply(t,phi),nl,gamma[mu])
            jw=I*sum((phi[a].conj()*t[a][bb]*dpcov[mu][bb]-dpcov[mu][a].conj()*t[a][bb]*phi[bb] for a,bb in product(range(2),repeat=2)),Q())
            literal-=g2*ji*jw/6
    cov_div=scalar_div+I*(b(kl,phi,nl)-b(ls,phi,an))
    check('full covariant fermion-current divergence jets',direct_div,cov_div)
    check('norm derivative includes both scalar slots',sum((j[mu]*dr[mu] for mu in range(4)),Q()),bop+dop)
    check('literal scalar gauge currents match three-tensor form',literal,I*(-2*g2*aop+(g2-gp2)*bop+(g2+gp2)*dop)/24)
    source=rr*(b(sb,phi,nl)-b(ls,phi,[add(msum((scale(m[i][j],nr[j]) for j in range(ns)),4,1),sr[i]) for i in range(ns)]))
    residual=rr*(b(eb,phi,nl)-b(ls,phi,er))
    totaldiv=sum((j[mu]*dr[mu] for mu in range(4)),Q())+rr*direct_div
    reduced=-I*(3*g2+gp2)*aop/24-I*gp2*bop/12+d*source
    check('complete scalar-current offshell IBP sourced residual',literal,reduced+I*d*totaldiv+d*residual)
    if trial==0:
        reject('opposite sterile EOM source sign',literal,reduced+I*d*totaldiv+d*residual+2*d*rr*b(ls,phi,sr))
        reject('drop differentiated fermion composite in IBP',literal,reduced+d*residual)
        reject('wrong lepton-source sign',literal,reduced+I*d*totaldiv+d*residual-2*d*rr*b(sb,phi,nl))
        wrong_scalar_div=sum((b(ls,dphi[mu],nl,gamma[mu]) for mu in range(4)),Q())
        reject('drop scalar connection in current divergence',direct_div,wrong_scalar_div+I*(b(kl,phi,nl)-b(ls,phi,an)))
        reject('discard scalar derivative remainder A',literal,reduced+I*d*totaldiv+d*residual+I*(3*g2+gp2)*aop/24)


def gv(a,v):return [sum((v[j]*a[i][j] for j in range(4)),Grass()) for i in range(4)]
def gb(v):return [sum((v[j].dagger()*gamma[0][j][i] for j in range(4)),Grass()) for i in range(4)]
def gc(v):return gv(C,gb(v))
def gadd(a,b):return [x+y for x,y in zip(a,b)]
def gscale(c,v):return [x*c for x in v]
def gsum(vs):
    z=[Grass() for s in range(4)]
    for v in vs:z=gadd(z,v)
    return z
def gdot(a,b):return sum((x*y for x,y in zip(a,b)),Grass())
def generators(offset):return [Grass({(offset+2*s,):1}) for s in range(4)]


na,ns,ne=2,3,2
phi=(Q(2,1),Q(-1,2));hc=(-phi[1],phi[0]);rr=sum((p.conj()*p for p in phi),Q())
f=rm(na,ns);y=rm(na,ns);ye=rm(na,ne);mass=mat(ns,ns,lambda i,j:(2,3,5)[i] if i==j else 0)
l=[[gv(PL,generators(16*al+8*a)) for a in range(2)] for al in range(na)]
barl=[[gb(l[al][a]) for a in range(2)] for al in range(na)]
nn=[gadd(generators(64+8*i),gc(generators(64+8*i))) for i in range(ns)]
nl=[gv(PL,n) for n in nn];nr=[gv(PR,n) for n in nn];barnr=[gb(n) for n in nr]
eright=[gv(PR,generators(128+8*i)) for i in range(ne)];bare=[gb(e) for e in eright]
xbar=[gsum(gscale(phi[a],barl[al][a]) for a in range(2)) for al in range(na)]
xc=[gv(C,xbar[al]) for al in range(na)]
sl=[gsum(gscale(y[al][i].conj()*phi[a].conj(),l[al][a]) for al in range(na) for a in range(2)) for i in range(ns)]
sr=[gc(s) for s in sl]
for i in range(ns):check('literal Majorana Yukawa source conjugation',sr[i],gsum(gscale(y[al][i],xc[al]) for al in range(na)))
sbar=[[gsum([gscale(y[al][j].conj()*phi[a].conj(),barnr[j]) for j in range(ns)]+[gscale(ye[al][r].conj()*hc[a],bare[r]) for r in range(ne)]) for a in range(2)] for al in range(na)]
raw=Grass();charged=Grass()
for al,i in product(range(na),range(ns)):
    row=gsum(gscale(phi[a],sbar[al][a]) for a in range(2))
    rcol=gsum([gscale(mass[i][j],nr[j]) for j in range(ns)]+[sr[i]])
    raw+=(gdot(row,nl[i])-gdot(xbar[al],rcol))*(rr*f[al][i])
    charged+=sum((gdot(bare[r],nl[i])*(rr*f[al][i]*ye[al][r].conj()*hc[a]*phi[a]) for r in range(ne) for a in range(2)),Grass())
hf=mul(dagger(y),f);k=mul(f,tuple(zip(*y)));fm=mul(f,mass)
sterile=sum((gdot(barnr[j],nl[i])*(rr*rr*hf[j][i]) for j,i in product(range(ns),repeat=2)),Grass())
mixed=sum((gdot(xbar[al],nr[i])*(rr*fm[al][i]) for al in range(na) for i in range(ns)),Grass())
active=sum((gdot(xbar[al],xc[be])*(rr*k[al][be]) for al,be in product(range(na),repeat=2)),Grass())
check('complete explicit scalar-current Yukawa and retained-mass source expansion',raw,sterile-mixed-active)
check('charged-lepton source vanishes after Hdagger phi contraction',charged,Grass())
for al,be in product(range(na),repeat=2):check('active scalar Grassmann flavor symmetry',gdot(xbar[al],xc[be]),gdot(xbar[be],xc[al]))
for j,i in product(range(ns),repeat=2):check('sterile Majorana scalar flavor symmetry',gdot(barnr[j],nl[i]),gdot(barnr[i],nl[j]))
check('active scalar source selects flavor symmetric part',active,sum((gdot(xbar[al],xc[be])*(rr*(k[al][be]+k[be][al])/2) for al,be in product(range(na),repeat=2)),Grass()))
check('sterile scalar source selects flavor symmetric part',sterile,sum((gdot(barnr[j],nl[i])*(rr*rr*(hf[j][i]+hf[i][j])/2) for j,i in product(range(ns),repeat=2)),Grass()))
active_pairs=sum((gdot(xbar[a],xc[a])*(rr*k[a][a]) for a in range(na)),Grass())+sum((gdot(xbar[a],xc[b])*(rr*(k[a][b]+k[b][a])) for a in range(na) for b in range(a+1,na)),Grass())
sterile_pairs=sum((gdot(barnr[i],nl[i])*(rr*rr*hf[i][i]) for i in range(ns)),Grass())+sum((gdot(barnr[j],nl[i])*(rr*rr*(hf[j][i]+hf[i][j])) for j in range(ns) for i in range(j+1,ns)),Grass())
check('active scalar independent flavor pair normalization',active,active_pairs)
check('sterile scalar independent flavor pair normalization',sterile,sterile_pairs)
cs=mat(na,na,lambda a,b:Q(a+b+1,a*b+2))
selectedf=mul(cs,tuple(tuple(z.conj() for z in row) for row in y));selectedh=mul(dagger(y),selectedf)
check('selected symmetric seed gives symmetric sterile scalar source',selectedh,tuple(zip(*selectedh)))
reject('active scalar incorrectly uses antisymmetric flavor projection',active,sum((gdot(xbar[al],xc[be])*(rr*(k[al][be]-k[be][al])/2) for al,be in product(range(na),repeat=2)),Grass()))
reject('drop paired active Majorana source',raw,sterile-mixed)
reject('opposite retained-mass source sign',raw,sterile+mixed-active)
reject('extra mixed Majorana half in source coefficient',raw,sterile-mixed*R(1,2)-active)


# Local field substitution implementing the new scalar residual. The
# Majorana variation is independently checked in the exterior algebra.
z=generators(192);en=gadd(z,gc(z));eR=gv(PR,en)
xx=gscale(-R(7,13),gsum(gscale(f[al][0].conj()*phi[a].conj(),l[al][a]) for al in range(na) for a in range(2)))
xreal=gadd(xx,gc(xx));packet=gdot(gb(xx),eR)
check('extra scalar map Majorana completed residual',packet+packet.dagger(),gdot(gb(xreal),en))
for a,b in zip(xreal,gc(xreal)):check('extra scalar map preserves Majorana reality',a,b)
reject('half the extra scalar Majorana field map',packet+packet.dagger(),gdot(gb(xreal),en)*R(1,2))


result={
    'status':'passed','python':sys.version.split()[0],
    'input_candidate':'5f9889b76661faa155c681a5ecfac5d036398bc8',
    'arithmetic':'Exact Gaussian rational weak polynomials, connection jets and finite exterior algebra; standard library only.',
    'positive_check_count':sum(counts.values()),'checks_by_category':dict(sorted(counts.items())),
    'negative_control_count':len(mutations),'negative_controls_rejected':sorted(mutations),
    'weak_covariance_systems':weak_receipts,
    'trial_configuration':{'generic_offshell_covariant_jet_trials':8,'active_flavors':2,'sterile_flavors':3,'charged_lepton_flavors':2,'complete_weak_polynomial_coefficients':True,'full_Grassmann_source_expansion':True},
    'scope':'Complete selected scalar gauge-current subpacket transport in a specified representative; not full derivative-EFT closure or a physical Wilson coefficient/beta function.'
}
Path(__file__).with_name('SCALAR_CURRENT_REDUCTION_RESULTS.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps(result,indent=2,sort_keys=True))
