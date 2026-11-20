#!/usr/bin/env python3
"""Zero-N source projection and scalar/derivative flavor statistics.

Standalone exact Gaussian rational Fourier/exterior algebra. Literal
neutrino sources enter every term of the existing sourced representative
before projection. No project imports or saved coefficient imports.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import time


class Q:
    __slots__=('r','i')
    def __init__(self,r=0,i=0):
        if isinstance(r,Q): self.r,self.i=r.r,r.i
        else: self.r,self.i=F(r),F(i)
    def __add__(self,o):
        if isinstance(o,G): return NotImplemented
        o=Q(o); return Q(self.r+o.r,self.i+o.i)
    __radd__=__add__
    def __neg__(self): return Q(-self.r,-self.i)
    def __sub__(self,o): return self+-Q(o)
    def __rsub__(self,o): return Q(o)+-self
    def __mul__(self,o):
        if isinstance(o,G): return NotImplemented
        o=Q(o); return Q(self.r*o.r-self.i*o.i,self.r*o.i+self.i*o.r)
    __rmul__=__mul__
    def __bool__(self): return bool(self.r or self.i)
    def __eq__(self,o):
        o=Q(o); return self.r==o.r and self.i==o.i
    def star(self): return Q(self.r,-self.i)


I=Q(0,1); ETA=(1,-1,-1,-1); ZERO=(0,0,0,0)


class G:
    """Finite Fourier exterior algebra, truncated at kappa^1 and fermion degree4."""
    def __init__(self,a=None):
        self.a={k:Q(v) for k,v in (a or {}).items() if Q(v) and k[1]<=1 and len(k[2])<=4}
    @staticmethod
    def coerce(x): return x if isinstance(x,G) else G({(ZERO,0,()):Q(x)})
    def __add__(self,o):
        o=G.coerce(o); out=dict(self.a)
        for k,v in o.a.items(): out[k]=out.get(k,Q())+v
        return G(out)
    __radd__=__add__
    def __neg__(self): return G({k:-v for k,v in self.a.items()})
    def __sub__(self,o): return self+-G.coerce(o)
    def __rsub__(self,o): return G.coerce(o)+-self
    def __mul__(self,o):
        o=G.coerce(o); out={}
        for (p,k,a),v in self.a.items():
            for (q,l,b),w in o.a.items():
                if k+l>1 or len(a)+len(b)>4 or set(a)&set(b): continue
                sign=(-1)**sum(x>y for x in a for y in b)
                key=(tuple(x+y for x,y in zip(p,q)),k+l,tuple(sorted(a+b)))
                out[key]=out.get(key,Q())+sign*v*w
        return G(out)
    __rmul__=__mul__
    def __eq__(self,o): return self.a==G.coerce(o).a
    def __bool__(self): return bool(self.a)
    def dagger(self):
        return G({(tuple(-x for x in p),k,a):v.star()*((-1)**(len(a)*(len(a)-1)//2)) for (p,k,a),v in self.a.items()})
    def deriv(self,mu): return G({(p,k,a):v*I*ETA[mu]*p[mu] for (p,k,a),v in self.a.items()})
    def shift(self): return G({(p,k+1,a):v for (p,k,a),v in self.a.items()})
    def order(self,k): return G({(p,0,a):v for (p,j,a),v in self.a.items() if j==k})
    def quadratic(self): return G({k:v for k,v in self.a.items() if len(k[2])==2})
    def integral(self): return G({k:v for k,v in self.a.items() if k[0]==ZERO})
    def without(self,ids): return G({k:v for k,v in self.a.items() if not (set(k[2])&ids)})


def theta(i,p): return G({(p,0,(i,)):1})
def vvadd(*vectors): return [sum((v[j] for v in vectors),G()) for j in range(4)]
def vvscale(c,v): return [c*x for x in v]
def mmv(m,v): return [sum((m[i][j]*v[j] for j in range(4)),G()) for i in range(4)]
def mm(a,b): return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(4)),Q()) for j in range(4)) for i in range(4))
def ms(c,a): return tuple(tuple(c*x for x in row) for row in a)
def mt(a): return tuple(tuple(a[j][i] for j in range(4)) for i in range(4))
def ma(a,b): return tuple(tuple(x+y for x,y in zip(ar,br)) for ar,br in zip(a,b))
def eye(n): return tuple(tuple(Q(i==j) for j in range(n)) for i in range(n))


PAULI=(((Q(),Q(1)),(Q(1),Q())),((Q(),-I),(I,Q())),((Q(1),Q()),(Q(),Q(-1))))
G0=tuple(tuple(Q((i<2 and j==i+2) or (i>=2 and j==i-2)) for j in range(4)) for i in range(4))
def spatial(t): return tuple(tuple(t[i][j-2] if i<2<=j else -t[i-2][j] if j<2<=i else Q() for j in range(4)) for i in range(4))
GAM=(G0,)+tuple(spatial(t) for t in PAULI)
CH=ms(I,mm(GAM[2],GAM[0])); BC=mm(CH,mt(G0))


def cc(v): return mmv(BC,[x.dagger() for x in v])
def bar(v):
    d=[x.dagger() for x in v]
    return [sum((d[i]*G0[i][j] for i in range(4)),G()) for j in range(4)]
def dot(row,col): return sum((a*b for a,b in zip(row,col)),G())
def A(v): return [sum((I*GAM[mu][j][k]*v[k].deriv(mu) for mu in range(4) for k in range(4)),G()) for j in range(4)]
def P0(v): return [sum((-ETA[mu]*x.deriv(mu).deriv(mu) for mu in range(4)),G()) for x in v]
def PP(v,m2): return vvadd(P0(v),vvscale(m2,v))
def Krow(row): return [sum((-I*row[j].deriv(mu)*GAM[mu][j][k] for mu in range(4) for j in range(4)),G()) for k in range(4)]
def hc(a): return a+a.dagger()
def integrate(a): return a.quadratic().integral()
def flmix(matrix,fields):
    return [vvadd(*(vvscale(c,fields[j]) for j,c in enumerate(row))) for row in matrix]
def adjmat(a): return tuple(tuple(a[j][i].star() for j in range(len(a))) for i in range(len(a[0])))


def mmrect(a,b):
    return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(len(b))),Q())
                       for j in range(len(b[0]))) for i in range(len(a)))
def mstar(a): return tuple(tuple(x.star() for x in row) for row in a)
def mtranspose(a): return tuple(tuple(a[j][i] for j in range(len(a))) for i in range(len(a[0])))
def msym(a): return tuple(tuple((a[i][j]+a[j][i])*F(1,2) for j in range(len(a))) for i in range(len(a)))
def mherm(a): return tuple(tuple((a[i][j]+a[j][i].star())*F(1,2) for j in range(len(a))) for i in range(len(a)))
def cfield(c,p=ZERO): return G({(p,0,()):c})
def covfield(v,mu,aa):
    return [v[a].deriv(mu)+I*sum((aa[mu][a][b]*v[b] for b in range(2)),G()) for a in range(2)]
def covrow(row,mu,aa):
    return [row[a].deriv(mu)-I*sum((row[b]*aa[mu][b][a] for b in range(2)),G()) for a in range(2)]


def make_fields(na,ns,p):
    nextgen=0
    def one():
        nonlocal nextgen
        v=[G() for _ in range(4)]
        for sign in (1,-1):
            momentum=tuple(sign*x for x in p)
            for spin in range(2):
                v[spin]+=theta(nextgen,momentum)+I*theta(nextgen+1,momentum)
                nextgen+=2
        return v
    ll=[[one() for w in range(2)] for a in range(na)]
    first_n=nextgen
    nn=[one() for i in range(ns)]
    nids=set(range(first_n,nextgen))
    return ll,nn,nids


def composites(ll,phi):
    xx=[vvadd(*(vvscale(phi[w].dagger(),l[w]) for w in range(2))) for l in ll]
    return xx,[bar(v) for v in xx],[cc(v) for v in xx]


def setup(na,p):
    ll,nl,nids=make_fields(na,2,p)
    phi=[cfield(Q(1,1))+cfield(Q(2,-1),(1,1,0,1))+cfield(Q(-1,3),(-1,-1,0,-1)),
         cfield(Q(-1,2))+cfield(Q(1,1),(-1,1,0,0))]
    yy=((Q(1,2),Q(2,-1)),(Q(3,1),Q(1,-1)),(Q(2,-2),Q(4,1)))[:na]
    seed=tuple(tuple(Q(i+j+1,(1 if i==j else -1)) for j in range(na)) for i in range(na))
    ff=mmrect(seed,mstar(yy))
    # Constant noncommuting Hermitian weak connections. Curvature is derived,
    # not copied from a previously simplified strength packet.
    aa=tuple(tuple(tuple(cfield(x*F(mu+1,5)) for x in row) for row in PAULI[mu%3]) for mu in range(4))
    ffcal=tuple(tuple(tuple(tuple(I*sum((aa[mu][a][c]*aa[nu][c][b]-aa[nu][a][c]*aa[mu][c][b]
                               for c in range(2)),G()) for b in range(2)) for a in range(2))
                         for nu in range(4)) for mu in range(4))
    return ll,nl,nids,phi,yy,ff,aa,ffcal


SIG=tuple(tuple(ms(I*F(1,2),ma(mm(GAM[mu],GAM[nu]),ms(-1,mm(GAM[nu],GAM[mu]))))
                for nu in range(4)) for mu in range(4))
def triple(mu,nu,rho):
    inds=(mu,nu,rho); out=ms(0,G0)
    for p in ((0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,0,1),(2,1,0)):
        sign=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
        out=ma(out,ms(F(sign,6),mm(mm(GAM[inds[p[0]]],GAM[inds[p[1]]]),GAM[inds[p[2]]])))
    return out
TRIPLE=tuple(tuple(tuple(triple(mu,nu,rho) for rho in range(4)) for nu in range(4)) for mu in range(4))


def literal_sources(ll,nl,phi,yy):
    nr=[cc(n) for n in nl]
    sl=[[vvadd(*(vvscale(yy[a][j]*phi[w],nr[j]) for j in range(len(nl)))) for w in range(2)] for a in range(len(ll))]
    sp=[-sum((yy[a][j].star()*dot(bar(nr[j]),ll[a][w]) for a,j in product(range(len(ll)),range(len(nl)))),G()) for w in range(2)]
    xx,_,_=composites(ll,phi)
    sln=flmix(adjmat(yy),xx); sr=[cc(v) for v in sln]
    return nr,sl,sp,sln,sr


def literal_representative(ll,nl,phi,yy,ff,aa,fc,masses,m2,lam,g2=F(7,5),b2=F(3,4),extra_sl=None,extra_sp=None):
    na,ns=len(ll),len(nl)
    nr,sl,sp,sln,sr=literal_sources(ll,nl,phi,yy)
    if extra_sl is not None:
        sl=[[vvadd(sl[a][w],extra_sl[a][w]) for w in range(2)] for a in range(na)]
    if extra_sp is not None: sp=[sp[w]+extra_sp[w] for w in range(2)]
    rr=[vvadd(vvscale(masses[i],nr[i]),sr[i]) for i in range(ns)]
    tt=[vvadd(vvscale(masses[i]**2,nl[i]),vvscale(masses[i],sln[i]),A(sr[i])) for i in range(ns)]
    uu=[vvadd(vvscale(masses[i]**3,nr[i]),vvscale(masses[i]**2,sr[i]),vvscale(masses[i],A(sln[i])),A(A(sr[i]))) for i in range(ns)]
    rho=sum((p.dagger()*p for p in phi),G()); c=2*m2*F(1,3)+lam*F(1,3)*rho
    terms={name:G() for name in ('U','cR','KbarS_R','Sphi_R','barS_T','m2barS_N','Vh','Eh','currents','strength_R')}
    dphi=[covfield(phi,mu,aa) for mu in range(4)]
    # K bar S is implemented as the full charged row derivative, then gamma.
    kbars=[]
    for a in range(na):
        brows=[bar(sl[a][w]) for w in range(2)]
        kr=[[G() for s in range(4)] for w in range(2)]
        for mu in range(4):
            for spin in range(4):
                derivative=covrow([brows[w][spin] for w in range(2)],mu,aa)
                for w,s in product(range(2),range(4)): kr[w][s]+=-I*derivative[w]*GAM[mu][spin][s]
        kbars.append(kr)
    for a,i,w in product(range(na),range(ns),range(2)):
        f=ff[a][i]; row=bar(ll[a][w])
        terms['U']+=f*F(5,6)*phi[w]*dot(row,uu[i])
        terms['cR']+=f*c*phi[w]*dot(row,rr[i])
        terms['KbarS_R']+=f*F(1,3)*phi[w]*dot(kbars[a][w],rr[i])
        terms['Sphi_R']-=f*F(1,6)*sp[w]*dot(row,rr[i])
        terms['barS_T']-=f*F(1,6)*phi[w]*dot(bar(sl[a][w]),tt[i])
        terms['m2barS_N']+=f*m2*F(1,2)*phi[w]*dot(bar(sl[a][w]),nl[i])
        for b,mu,nu in product(range(2),range(4),range(4)):
            terms['strength_R']-=f*F(1,12)*fc[mu][nu][w][b]*phi[b]*dot(row,mmv(SIG[mu][nu],rr[i]))
            terms['Vh']+=f*F(1,2)*fc[mu][nu][w][b]*ETA[nu]*dphi[nu][b]*dot(row,mmv(GAM[mu],nl[i]))
            for r in range(4): terms['Eh']+=f*F(1,4)*fc[mu][nu][w][b]*dphi[r][b]*dot(row,mmv(TRIPLE[mu][nu][r],nl[i]))
    # The actual retained scalar + left-lepton currents, with their complete
    # connection terms, have no N. Their product with j always retains N.
    tgens=tuple(tuple(tuple(x*F(1,2) for x in row) for row in t) for t in PAULI)
    for mu in range(4):
        j=sum((ff[a][i]*phi[w]*ETA[mu]*dot(bar(ll[a][w]),mmv(GAM[mu],nl[i]))
               for a,i,w in product(range(na),range(ns),range(2))),G())
        jy=-F(1,2)*sum((dot(bar(ll[a][w]),mmv(GAM[mu],ll[a][w])) for a,w in product(range(na),range(2))),G())
        jy+=-I*F(1,2)*ETA[mu]*sum((phi[w].dagger()*dphi[mu][w]-dphi[mu][w].dagger()*phi[w] for w in range(2)),G())
        terms['currents']+=b2*F(1,12)*j*jy
        for t in tgens:
            ji=sum((ff[a][i]*t[w][b]*phi[b]*ETA[mu]*dot(bar(ll[a][w]),mmv(GAM[mu],nl[i]))
                    for a,i,w,b in product(range(na),range(ns),range(2),range(2))),G())
            jf=sum((t[w][b]*dot(bar(ll[a][w]),mmv(GAM[mu],ll[a][b])) for a,w,b in product(range(na),range(2),range(2))),G())
            jf+=I*ETA[mu]*sum((phi[w].dagger()*t[w][b]*dphi[mu][b]-dphi[mu][w].dagger()*t[w][b]*phi[b] for w,b in product(range(2),repeat=2)),G())
            terms['currents']-=g2*F(1,6)*ji*jf
    return terms


def scalar_packet(xx,bx,xc,matrix,operation='scalar',rho=None):
    out=G()
    for a,b in product(range(len(xx)),repeat=2):
        right=P0(xc[b]) if operation=='box' else A(xx[b]) if operation=='dirac' else xc[b]
        out+=matrix[a][b]*dot(bx[a],right)
    return (rho*out if rho is not None else out)


def strength_packet(ll,phi,xc,kk,fc):
    return sum((kk[a][b]*fc[mu][nu][w][v]*phi[v]*dot(bar(ll[a][w]),mmv(SIG[mu][nu],xc[b]))
                for a,b,w,v,mu,nu in product(range(len(ll)),range(len(ll)),range(2),range(2),range(4),range(4))),G())


def compact_projection(ll,phi,yy,ff,fc,masses,m2,lam):
    xx,bx,xc=composites(ll,phi)
    kk=mmrect(ff,mtranspose(yy))
    km2=mmrect(tuple(tuple(ff[a][i]*masses[i]**2 for i in range(len(masses))) for a in range(len(ll))),mtranspose(yy))
    gg=mmrect(tuple(tuple(ff[a][i]*masses[i] for i in range(len(masses))) for a in range(len(ll))),adjmat(yy))
    rho=sum((p.dagger()*p for p in phi),G())
    pieces={'mass':F(5,6)*scalar_packet(xx,bx,xc,km2),
            'box':F(5,6)*scalar_packet(xx,bx,xc,kk,'box'),
            'dirac':F(5,6)*scalar_packet(xx,bx,xc,gg,'dirac'),
            'm2':2*m2*F(1,3)*scalar_packet(xx,bx,xc,kk),
            'rho':lam*F(1,3)*scalar_packet(xx,bx,xc,kk,rho=rho),
            'strength':-F(1,12)*strength_packet(ll,phi,xc,kk,fc)}
    return pieces,(xx,bx,xc,kk,km2,gg,rho)


def source_map_active_variation(ll,phi,yy,ff,fc,masses,m2,lam,aa=None,extra_sl=None,extra_sp=None):
    """Substitute actual source-map sterile shift into reference Yukawa at Nnew=0."""
    xx,bx,xc=composites(ll,phi)
    ux=flmix(adjmat(ff),xx); rho=sum((p.dagger()*p for p in phi),G())
    vnr=[]
    for i in range(len(masses)):
        x=vvadd(vvscale(F(5,6),P0(ux[i])),vvscale(F(5,6)*masses[i]**2,ux[i]),vvscale(2*m2*F(1,3)+lam*F(1,3)*rho,ux[i]))
        for a,w,b,mu,nu in product(range(len(ll)),range(2),range(2),range(4),range(4)):
            x=vvadd(x,vvscale(-F(1,12)*ff[a][i].star()*phi[b].dagger()*fc[mu][nu][b][w],mmv(SIG[mu][nu],ll[a][w])))
        y=vvscale(F(5,6)*masses[i],A(ux[i]))
        if extra_sl is not None:
            dphi=[covfield(phi,mu,aa) for mu in range(4)]
            for a,w,mu in product(range(len(ll)),range(2),range(4)):
                dsl=[covfield([extra_sl[a][v][spin] for v in range(2)],mu,aa)[w] for spin in range(4)]
                x=vvadd(x,vvscale(I*F(1,6)*ff[a][i].star()*phi[w].dagger(),mmv(GAM[mu],dsl)),
                          vvscale(-I*F(1,6)*ff[a][i].star()*dphi[mu][w].dagger(),mmv(GAM[mu],extra_sl[a][w])))
            for a,w in product(range(len(ll)),range(2)):
                y=vvadd(y,vvscale(-F(1,6)*ff[a][i].star()*masses[i]*phi[w].dagger(),extra_sl[a][w]))
        if extra_sp is not None:
            for a,w in product(range(len(ll)),range(2)):
                x=vvadd(x,vvscale(-F(1,6)*ff[a][i].star()*extra_sp[w].dagger(),ll[a][w]))
        vnr.append(vvscale(-1,vvadd(y,cc(x))))
    # Reference Yukawa -barX Y N_R+h.c.; no original N remains.
    return hc(-sum((yy[a][i]*dot(bx[a],vnr[i]) for a,i in product(range(len(ll)),range(len(masses)))),G()))


def original_scalar_current(ll,nl,phi,ff,aa):
    """Actual O1,O2,O3 combination with g^2=4, gprime^2=3."""
    rho=sum((v.dagger()*v for v in phi),G()); o1=o2=o3=G()
    for mu in range(4):
        dp=covfield(phi,mu,aa)
        brow=[sum((ff[a][i]*dot(bar(ll[a][w]),mmv(GAM[mu],nl[i])) for a,i in product(range(len(ll)),range(len(nl)))),G()) for w in range(2)]
        j=sum((brow[w]*phi[w] for w in range(2)),G())
        o1+=ETA[mu]*rho*sum((brow[w]*dp[w] for w in range(2)),G())
        o2+=ETA[mu]*j*sum((phi[w].dagger()*dp[w] for w in range(2)),G())
        o3+=ETA[mu]*j*sum((dp[w].dagger()*phi[w] for w in range(2)),G())
    return I*F(1,24)*(-8*o1+o2+7*o3)


def spectator_sources(ll,nl,phi,p):
    newll,newnl,_=make_fields(2,2,p)
    offset=1+max(g for group in ll+[nl] for v in group for x in v for key in x.a for g in key[2])
    def move(v): return [G({(mom,k,tuple(g+offset for g in word)):c for (mom,k,word),c in x.a.items()}) for x in v]
    er=cc(move(newll[0][0])); dr=cc(move(newll[0][1])); ur=cc(move(newnl[0]))
    ql=[move(v) for v in newll[1]]
    hh=[-phi[1].dagger(),phi[0].dagger()]
    ye=[Q(a+1,2-a) for a in range(len(ll))]; yd=Q(2,-1); yu=Q(1,2)
    esl=[[vvscale(ye[a]*hh[w],er) for w in range(2)] for a in range(len(ll))]
    eps=((0,1),(-1,0))
    esp=[]
    for b in range(2):
        value=sum((ye[a]*eps[w][b]*dot(bar(ll[a][w]),er) for a,w in product(range(len(ll)),range(2))),G())
        value+=sum((yd*eps[w][b]*dot(bar(ql[w]),dr) for w in range(2)),G())
        value-=yu.star()*dot(bar(ur),ql[b])
        esp.append(value)
    return esl,esp,hh


def spectator_remainder(ll,phi,yy,ff,aa,esl,esp):
    _,_,_,_,sr=literal_sources(ll,[[G() for s in range(4)] for j in range(len(yy[0]))],phi,yy)
    out=G()
    for a,i,w in product(range(len(ll)),range(len(sr)),range(2)):
        row=bar(esl[a][w]); krow=[G() for s in range(4)]
        for mu,spin in product(range(4),range(4)):
            derivative=covrow([bar(esl[a][v])[spin] for v in range(2)],mu,aa)[w]
            for s in range(4): krow[s]+=-I*derivative*GAM[mu][spin][s]
        out+=F(1,3)*ff[a][i]*phi[w]*dot(krow,sr[i])
        out-=F(1,6)*ff[a][i]*esp[w]*dot(bar(ll[a][w]),sr[i])
    return out


def rank(vectors):
    # Gaussian rational division is implemented locally to keep Q minimal.
    def divide(x,y):
        d=y.r*y.r+y.i*y.i
        return Q((x.r*y.r+x.i*y.i)/d,(x.i*y.r-x.r*y.i)/d)
    pivots={}
    for vector in vectors:
        v=dict(vector.a)
        while v:
            p=min(v)
            if p not in pivots:
                pivots[p]={k:divide(x,v[p]) for k,x in v.items()}; break
            lead=v[p]
            for k,c in pivots[p].items(): v[k]=v.get(k,Q())-lead*c
            v={k:c for k,c in v.items() if c}
    return len(pivots)


checks={}; mutations={}; degree_records=[]; rank_records=[]
def check(cat,label,a,b):
    if a!=b: raise RuntimeError('FAILED '+cat+': '+label)
    checks[cat]=checks.get(cat,0)+1
def reject(label,a,b):
    if a==b: raise RuntimeError('MUTATION SURVIVED '+label)
    mutations[label]=mutations.get(label,0)+1


def run():
    started=time.perf_counter(); masses=(F(2),F(5)); m2=F(7,3); lam=F(2,5)
    for na,p in ((1,(2,1,0,-1)),(2,(3,1,-1,1))):
        print('Active projection: %s active flavors' % na,flush=True)
        ll,nl,nids,phi,yy,ff,aa,fc=setup(na,p)
        terms=literal_representative(ll,nl,phi,yy,ff,aa,fc,masses,m2,lam)
        expected,inputs=compact_projection(ll,phi,yy,ff,fc,masses,m2,lam)
        direct=sum(terms.values(),G()).without(nids); compact=sum(expected.values(),G())
        check('literal_full_source_projection','all retained-N content removed',direct,compact)
        for name,term in terms.items():
            n_degrees=sorted(set(len(set(key[2])&nids) for key in term.a))
            degree_records.append({'active_flavors':na,'term':name,'retained_N_degrees':n_degrees,
                                   'fermion_degrees':sorted(set(len(key[2]) for key in term.a)),
                                   'zero_N_monomials':len(term.without(nids).a)})
            if name not in ('U','cR','strength_R'):
                check('vanishing_source_classes',name,term.without(nids),G())
        xx,bx,xc,kk,km2,gg,rho=inputs
        for name,mat,op,weight in (('scalar',kk,'scalar',None),('box',kk,'box',None),('mass',km2,'scalar',None),('rho',kk,'scalar',rho)):
            raw=scalar_packet(xx,bx,xc,mat,op,weight)
            sy=scalar_packet(xx,bx,xc,msym(mat),op,weight)
            check('scalar_flavor_statistics',name,raw.integral(),sy.integral())
        rawd=scalar_packet(xx,bx,xc,gg,'dirac')
        hermd=scalar_packet(xx,bx,xc,mherm(gg),'dirac')
        check('LNC_hermitian_transport','one Hermitian action, no repeated h.c.',hc(rawd).integral(),(2*hermd).integral())
        check('original_coordinate_source_variation','actual zero-N source-map Yukawa variation',source_map_active_variation(ll,phi,yy,ff,fc,masses,m2,lam).integral(),hc(compact).integral())
        # Map the Box source into the literal seed convention +1/2 C7 O7.
        delta_c7=tuple(tuple(F(5,3)*x for x in row) for row in msym(kk))
        check('literal_action_normalization','seed half factor',(hc(expected['box'])).integral(),(F(1,2)*hc(scalar_packet(xx,bx,xc,delta_c7,'box'))).integral())
        reject('spurious_extra_seed_half',hc(expected['box']).integral(),(F(1,4)*hc(scalar_packet(xx,bx,xc,delta_c7,'box'))).integral())
        reject('delete_retained_mass_source',hc(compact).integral(),hc(compact-expected['mass']-expected['dirac']).integral())
        reject('discard_active_projection_as_no_sterile_external_field',hc(compact).integral(),G())
        reject('derivative_only_on_lepton_not_phi',expected['box'].integral(),
               (F(5,6)*sum((kk[a][b]*phi[v]*dot(bx[a],P0(cc(ll[b][v])))
                            for a,b,v in product(range(na),range(na),range(2))),G())).integral())
        reject('omit_LNC_hermitian_conjugate',hc(rawd).integral(),rawd.integral())
        # Additional scalar-current map: project the old scalar subpacket vs
        # pull back the reference Yukawa with Nold=kappa d rho Fdagger X.
        d=F(7,24)  # a generic positive (g^2+gprime^2)/24 benchmark
        vnl=flmix(adjmat(ff),xx); vnl=[vvscale(d*rho,v) for v in vnl]
        vnr=[cc(v) for v in vnl]
        delta_y=hc(-sum((yy[a][i]*dot(bx[a],vnr[i]) for a,i in product(range(na),range(2))),G()))
        induced=-d*hc(scalar_packet(xx,bx,xc,kk,rho=rho))
        check('scalar_map_projection_noncommutation','actual shifted Yukawa',delta_y,induced)
        scalar_before=original_scalar_current(ll,nl,phi,ff,aa)
        check('scalar_map_projection_noncommutation','before-map scalar subpacket has N',scalar_before.without(nids),G())
        reject('scalar_subpacket_accidentally_zero_before_projection',scalar_before,G())
        reject('projection_commutes_with_scalar_field_map',delta_y,G())
        reject('wrong_scalar_map_active_sign',delta_y,-induced)
        if na>1:
            anti=tuple(tuple((kk[a][b]-kk[b][a])*F(1,2) for b in range(na)) for a in range(na))
            check('scalar_flavor_statistics','antisymmetric scalar zero',scalar_packet(xx,bx,xc,anti).integral(),G())
            check('scalar_flavor_statistics','antisymmetric Box zero after IBP',scalar_packet(xx,bx,xc,anti,'box').integral(),G())
            reject('antisymmetric_instead_of_symmetric_Box',expected['box'].integral(),scalar_packet(xx,bx,xc,anti,'box').integral())
            reject('weak_strength_reduced_by_flavor_symmetry_alone',expected['strength'], -F(1,12)*strength_packet(ll,phi,xc,msym(kk),fc))
        if na==1:
            # Explicit e,d,u spectator Yukawas: the main nu-only exhaustion
            # acquires these extra terms. N=0 never deletes these sources.
            esl,esp,hh=spectator_sources(ll,nl,phi,p)
            for a in range(na):
                check('SM_spectator_remainder','undifferentiated barS_e phi zero',
                      vvadd(*(vvscale(phi[w],bar(esl[a][w])) for w in range(2))),[G() for s in range(4)])
            nlzero=[[G() for s in range(4)] for i in range(2)]
            smterms=literal_representative(ll,nlzero,phi,yy,ff,aa,fc,masses,m2,lam,extra_sl=esl,extra_sp=esp)
            remainder=spectator_remainder(ll,phi,yy,ff,aa,esl,esp)
            check('SM_spectator_remainder','literal source additions after projection',sum(smterms.values(),G()),compact+remainder)
            check('SM_spectator_remainder','actual sterile source-map variation with spectators',
                  source_map_active_variation(ll,phi,yy,ff,fc,masses,m2,lam,aa,esl,esp).integral(),hc(compact+remainder).integral())
            reject('delete_all_SM_sources_when_N_is_zero',hc(compact+remainder).integral(),hc(compact).integral())
    # Independent operator rank controls in these algebraic field contents.
    for na in (1,2,3):
        ll,nl,nids,phi,yy,ff,aa,fc=setup(na,(3,1,0,-1))
        xx,bx,xc=composites(ll,phi)
        scalars=[]; boxes=[]; herm=[]
        for a in range(na):
            for b in range(a,na):
                mat=tuple(tuple(Q(i==a and j==b) for j in range(na)) for i in range(na))
                scalars.append(scalar_packet(xx,bx,xc,mat).integral())
                boxes.append(scalar_packet(xx,bx,xc,mat,'box').integral())
                mat_real=tuple(tuple(Q((i==a and j==b) or (i==b and j==a)) for j in range(na)) for i in range(na))
                herm.append(scalar_packet(xx,bx,xc,mat_real,'dirac').integral())
                if a!=b:
                    mat_imag=tuple(tuple(I*((i==a and j==b)-(i==b and j==a)) for j in range(na)) for i in range(na))
                    herm.append(scalar_packet(xx,bx,xc,mat_imag,'dirac').integral())
        check('component_operator_ranks','symmetric scalar',rank(scalars),na*(na+1)//2)
        check('component_operator_ranks','symmetric integrated Box',rank(boxes),na*(na+1)//2)
        check('component_operator_ranks','Hermitian LNC Dirac',rank(herm),na*na)
        rank_records.append({'active_flavors':na,'scalar_complex_rank':rank(scalars),'Box_complex_rank':rank(boxes),'Hermitian_LNC_real_rank':rank(herm)})
    result={'status':'PASS','positive_checks':sum(checks.values()),'checks_by_category':checks,
            'negative_checks':sum(mutations.values()),'mutation_rejections':mutations,
            'total_checks':sum(checks.values())+sum(mutations.values()),'term_N_degree_records':degree_records,
            'operator_ranks':rank_records,'arithmetic':'exact Gaussian rational Fourier exterior algebra, standard library only',
            'projection':'set retained-N component generators to zero; no sterile equation or integration used',
            'scope':['complete zero-N field-content projection of the declared literal-neutrino-sourced representative',
                     'scalar/Box flavor symmetry, LNC integrated Hermiticity, literal action half normalization',
                     'independent direct reference-Yukawa variation under the actual source map at Nnew=0',
                     'additional scalar-current map generates active rho scalar source after projection',
                     'explicit e,d,u spectator-Yukawa remainder checked separately; no general extra-interaction exhaustion',
                     'real positive diagonal retained M; complex-mass extension is separate',
                     'coefficient projection only; no physical beta, matching, finite amplitude, or new loop']}
    Path(__file__).with_name('ACTIVE_SOURCE_PROJECTION_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Exact active-source projection replay in %.3f seconds' % (time.perf_counter()-started),flush=True)
    print(json.dumps(result,indent=2))


if __name__=='__main__': run()


