#!/usr/bin/env python3
"""Actual neutral action pullbacks, nonlinear Hessian term, and map bridge.

Standalone exact Gaussian rational Fourier/exterior algebra. All Majorana
slots, Yukawa sources, source probes, and scalar potential are substituted
before extracting the quadratic, first-order coefficient. No project imports.
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
    """Finite Fourier exterior algebra, truncated at kappa^1 and degree2."""
    def __init__(self,a=None):
        self.a={k:Q(v) for k,v in (a or {}).items() if Q(v) and k[1]<=1 and len(k[2])<=2}
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
                if k+l>1 or len(a)+len(b)>2 or set(a)&set(b): continue
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


def sources(nu,nl,u,y):
    nr=[cc(v) for v in nl]
    sl=flmix(y,[vvscale(u,v) for v in nr])
    sln=flmix(adjmat(y),[vvscale(u,v) for v in nu])
    sr=[cc(v) for v in sln]
    return nr,sl,sln,sr


def action0(nu,nl,phi,phi_d,m2,lam,y,masses):
    total=G()
    for v in nu: total+=dot(bar(v),A(v))
    nr=[cc(v) for v in nl]
    for j in range(len(nl)):
        n=vvadd(nl[j],nr[j])
        total+=F(1,2)*dot(bar(n),vvadd(A(n),vvscale(-masses[j],n)))
    for a,j in product(range(len(nu)),range(len(nl))):
        total-=y[a][j]*phi*dot(bar(nu[a]),nr[j])
        total-=y[a][j].star()*phi_d*dot(bar(nr[j]),nu[a])
    rho=phi_d*phi
    total-=m2*rho+lam*rho*rho
    total+=sum((ETA[mu]*phi_d.deriv(mu)*phi.deriv(mu) for mu in range(4)),G())
    return total


def packet_original(nu,nl,u,m2,ff,collapsed=False):
    total=G()
    for a,j in product(range(len(nu)),range(len(nl))):
        row=bar(nu[a])
        if collapsed: term=dot(row,A(PP(nl[j],m2)))
        else:
            term=F(5,6)*dot(row,A(P0(nl[j])))
            term+=F(1,3)*dot(P0(row),A(nl[j]))
            term-=F(1,6)*dot(Krow(row),P0(nl[j]))
            term+=m2*F(1,2)*dot(row,A(nl[j]))
            term+=m2*F(1,2)*dot(Krow(row),nl[j])
        total+=u*ff[a][j]*term
    return hc(total)


def packet_reduced(nu,nl,u,m2,lam,ff,y,masses):
    nr,sl,sln,sr=sources(nu,nl,u,y)
    rr=[vvadd(vvscale(masses[j],nr[j]),sr[j]) for j in range(len(nl))]
    tt=[vvadd(vvscale(masses[j]**2,nl[j]),vvscale(masses[j],sln[j]),A(sr[j])) for j in range(len(nl))]
    uu=[vvadd(vvscale(masses[j]**3,nr[j]),vvscale(masses[j]**2,sr[j]),vvscale(masses[j],A(sln[j])),P0(sr[j])) for j in range(len(nl))]
    cv=2*m2*F(1,3)+lam*u*u*F(1,3)
    total=G()
    for a,j in product(range(len(nu)),range(len(nl))):
        term=F(5,6)*dot(bar(nu[a]),uu[j])+cv*dot(bar(nu[a]),rr[j])
        term+=F(1,3)*dot(Krow(bar(sl[a])),rr[j])
        term-=F(1,6)*dot(bar(sl[a]),tt[j])
        term+=m2*F(1,2)*dot(bar(sl[a]),nl[j])
        total+=u*ff[a][j]*term
    return hc(total)


def shifts(nu,nl,u,m2,lam,ff,y,masses,kind,sterile_factor=1,drop_yc=False,wrong_mass_order=False):
    na,ns=len(nu),len(nl)
    if kind=='C':
        dnu=[vvscale(-u*F(1,2),vvadd(*(vvscale(ff[a][j],PP(nl[j],m2)) for j in range(ns)))) for a in range(na)]
        dnl=[vvscale(-u*F(1,2),PP(vvadd(*(vvscale(ff[a][j].star(),nu[a]) for a in range(na))),m2)) for j in range(ns)]
        return dnu,dnl,G(),G()
    nr=[cc(v) for v in nl]; cv=2*m2*F(1,3)+lam*u*u*F(1,3)
    dnu=[vvscale(-u,vvadd(*(vvscale(ff[a][j],vvadd(vvscale(F(1,6),P0(nl[j])),vvscale(m2*F(1,2),nl[j]))) for j in range(ns)))) for a in range(na)]
    xs=[]; ys=[]
    for j in range(ns):
        x=vvadd(*(vvscale(u*ff[a][j].star(),vvadd(vvscale(F(5,6),P0(nu[a])),vvscale(cv+F(5,6)*masses[j]**2,nu[a]))) for a in range(na)))
        yr=vvadd(*(vvscale(F(5,6)*u*ff[a][j].star()*masses[a if wrong_mass_order else j],A(nu[a])) for a in range(na)))
        for a,k in product(range(na),range(ns)):
            x=vvadd(x,vvscale(u*u*F(1,6)*ff[a][j].star()*y[a][k],A(nr[k])))
            yr=vvadd(yr,vvscale(-u*u*F(1,6)*ff[a][j].star()*masses[a if wrong_mass_order else j]*y[a][k],nr[k]))
        xs.append(x); ys.append(yr)
    dnl=[vvscale(-sterile_factor,vvadd(xs[j],([G()]*4 if drop_yc else cc(ys[j])))) for j in range(ns)]
    t=sum((ff[a][j]*F(1,6)*dot(bar(nu[a]),A(nl[j])) for a,j in product(range(na),range(ns))),G())
    return dnu,dnl,-t.dagger(),-t


def add_shifts(a,b,sign=1):
    return ([vvadd(x,vvscale(sign,y)) for x,y in zip(a[0],b[0])],
            [vvadd(x,vvscale(sign,y)) for x,y in zip(a[1],b[1])],a[2]+sign*b[2],a[3]+sign*b[3])
def pull0(nu,nl,u,m2,lam,y,masses,v,scalar=True):
    nu1=[vvadd(f,[x.shift() for x in df]) for f,df in zip(nu,v[0])]
    nl1=[vvadd(f,[x.shift() for x in df]) for f,df in zip(nl,v[1])]
    phi=G.coerce(u)+(v[2].shift() if scalar else G())
    phi_d=G.coerce(u)+(v[3].shift() if scalar else G())
    return action0(nu1,nl1,phi,phi_d,m2,lam,y,masses)


def make_fields(na,ns,p):
    registry=[]; nextgen=0
    def one(kind,flavor):
        nonlocal nextgen
        v=[G() for _ in range(4)]
        for k,sgn in enumerate((1,-1)):
            momentum=tuple(sgn*x for x in p)
            for spin in range(2):
                re,im=nextgen,nextgen+1; nextgen+=2
                v[spin]+=theta(re,momentum)+I*theta(im,momentum)
                registry.append((kind,flavor,spin,momentum,re,im))
        return v
    nu=[one('L',a) for a in range(na)]; nl=[one('N',j) for j in range(ns)]
    return nu,nl,registry,nextgen


def jacobian(v,registry,n):
    out=[[Q() for j in range(n)] for i in range(n)]
    for kind,f,s,p,re,im in registry:
        field=v[0 if kind=='L' else 1][f][s]
        for (mom,k,word),c in field.a.items():
            if mom!=p or k or len(word)!=1: continue
            out[re][word[0]]=Q(c.r); out[im][word[0]]=Q(c.i)
    return out
def hessian(poly,n):
    out=[[Q() for j in range(n)] for i in range(n)]
    for (p,k,word),c in integrate(poly).a.items():
        if k or len(word)!=2: raise ValueError('Hessian extraction order')
        i,j=word; out[i][j]+=c; out[j][i]-=c
    return out
def congruence_variation(h,jac):
    n=len(h)
    return [[sum((jac[k][i]*h[k][j]+h[i][k]*jac[k][j] for k in range(n)),Q()) for j in range(n)] for i in range(n)]
def hadd(a,b): return [[x+y for x,y in zip(ar,br)] for ar,br in zip(a,b)]


checks={}; mutations={}; trials=[]
def check(cat,label,a,b):
    if a!=b: raise RuntimeError('FAILED '+cat+': '+label)
    checks[cat]=checks.get(cat,0)+1
def reject(label,a,b):
    if a==b: raise RuntimeError('MUTATION SURVIVED '+label)
    mutations[label]=mutations.get(label,0)+1


def run():
    started=time.perf_counter()
    ff=((Q(1,2),Q(2,-1)),(Q(-1,1),Q(3,2)))
    yy=((Q(2,-1),Q(1,3)),(Q(-2,2),Q(4,-1)))
    masses=(F(2),F(5)); u=F(3,2); lam=F(2,5)
    for p,m2 in (((2,1,-1,1),F(7,3)),((3,2,1,-1),-F(4,3)),((2,0,1,1),-2*lam*u*u)):
        nu,nl,reg,n=make_fields(2,2,p)
        base=integrate(action0(nu,nl,G.coerce(u),G.coerce(u),m2,lam,yy,masses))
        original=integrate(packet_original(nu,nl,u,m2,ff))
        reduced=integrate(packet_reduced(nu,nl,u,m2,lam,ff,yy,masses))
        check('literal_packet_projection','five nonvanishing of six placements',original,integrate(packet_original(nu,nl,u,m2,ff,collapsed=True)))
        vs=shifts(nu,nl,u,m2,lam,ff,yy,masses,'S')
        vc=shifts(nu,nl,u,m2,lam,ff,yy,masses,'C')
        bridge=add_shifts(vc,vs,-1)
        pulled_s=integrate(pull0(nu,nl,u,m2,lam,yy,masses,vs).order(1))+original
        pulled_c=integrate(pull0(nu,nl,u,m2,lam,yy,masses,vc).order(1))+original
        pulled_b=integrate(pull0(nu,nl,u,m2,lam,yy,masses,bridge).order(1))+reduced
        check('direct_actual_action_substitution','source representative versus original pullback',pulled_s,reduced)
        check('direct_actual_action_substitution','source-to-normalized coordinate bridge',pulled_b,pulled_c)
        # Scalar potential substituted independently, not appended as a correction.
        td=-(u*(m2+2*lam*u*u))*hc(vs[3])
        tad=integrate(td)
        only_j_s=integrate(pull0(nu,nl,u,m2,lam,yy,masses,vs,scalar=False).order(1))
        check('nonlinear_scalar_hessian','potential substitution supplies tadpole',pulled_s-original,only_j_s+tad)
        js=jacobian(vs,reg,n); jc=jacobian(vc,reg,n); jb=jacobian(bridge,reg,n)
        h0=hessian(base,n); hp=hessian(original,n); ht=hessian(tad,n)
        hs=hessian(reduced,n); hc0=hessian(pulled_c,n)
        check('explicit_hessian_congruence','source Hessian including gradient term',hs,hadd(hadd(hp,congruence_variation(h0,js)),ht))
        check('explicit_hessian_congruence','normalization Hessian',hc0,hadd(hp,congruence_variation(h0,jc)))
        check('explicit_hessian_congruence','bridge Hessian including opposite gradient',hc0,hadd(hadd(hs,congruence_variation(h0,jb)),[[-c for c in row] for row in ht]))
        # Check actual Majorana/chiral completion and derivative charge conjugation.
        for j in range(2):
            check('graded_majorana','charge-conjugation square',cc(cc(nl[j])),nl[j])
            check('graded_majorana','A commutes with charge conjugation',cc(A(nl[j])),A(cc(nl[j])))
            check('graded_majorana','constrained full sterile shift',cc(vvadd(vs[1][j],cc(vs[1][j]))),vvadd(vs[1][j],cc(vs[1][j])))
            check('graded_majorana','integrated row derivative identity',integrate(dot(Krow(bar(nu[j])),nl[j])),integrate(dot(bar(nu[j]),A(nl[j]))))
        # Normalization cancels the selected kinetic packet, retaining mass terms.
        def kinetic(nu0,nl0):
            return sum((dot(bar(v),A(v)) for v in nu0),G())+sum((F(1,2)*dot(bar(vvadd(v,cc(v))),A(vvadd(v,cc(v)))) for v in nl0),G())
        nuc=[vvadd(a,[x.shift() for x in b]) for a,b in zip(nu,vc[0])]
        nlc=[vvadd(a,[x.shift() for x in b]) for a,b in zip(nl,vc[1])]
        check('normalization_transport','selected kinetic cancellation',integrate(kinetic(nuc,nlc).order(1))+original,G())
        reject('omit_induced_mass_terms',pulled_c,G())
        # Probe rows use independent generators and both Fourier modes.
        probe_nu,probe_nl,_,_=make_fields(2,2,p)
        def offset(vv):
            return [G({(mom,k,tuple(g+n for g in word)):c for (mom,k,word),c in x.a.items()}) for x in vv]
        probe_nu=[offset(v) for v in probe_nu]; probe_nl=[offset(v) for v in probe_nl]
        def probe(v):
            # Rows have the chirality required to pair a left column; their
            # components are independent external odd variables, not S_L.
            return integrate(hc(sum((dot(bar(cc(probe_nu[a])),v[0][a]) for a in range(2)),G())+
                                sum((dot(bar(cc(probe_nl[j])),v[1][j]) for j in range(2)),G())))
        check('independent_source_transport','source bridge composition',probe(add_shifts(bridge,vs)),probe(vc))
        reject('omit_source_transport',probe(vc),G())
        # Even external sources are also pulled back. At this neutral
        # background only j3=j/2 is nonzero; charged-current probes belong to
        # the full interacting bridge in the derivation, not this Hessian.
        gp,g=F(3,5),F(4,7)
        jmu=[sum((u*ff[a][j]*dot(bar(nu[a]),vvscale(ETA[mu],mmv(GAM[mu],nl[j])))
                   for a,j in product(range(2),range(2))),G()) for mu in range(4)]
        source_b=integrate(sum((F(mu+1,3)*gp*(-F(1,2))*F(1,6)*hc(jmu[mu]) for mu in range(4)),G()))
        source_w=integrate(sum((F(2*mu+1,5)*g*F(1,12)*hc(jmu[mu]) for mu in range(4)),G()))
        source_phi=integrate(F(5,7)*(vs[2]+vs[3]))
        source_all=source_b+source_w+source_phi
        reject('omit_gauge_composite_source',source_all,source_phi)
        reject('omit_scalar_composite_source',source_all,source_b+source_w)
        # A separately declared fixed linear scalar source can support the
        # chosen off-shell background, but its pullback cancels the tadpole.
        supporting_action=u*(m2+2*lam*u*u)*(G.coerce(u)+vs[2].shift()+G.coerce(u)+vs[3].shift())
        supporting=integrate(supporting_action.order(1))
        check('independent_source_transport','supporting scalar source is transformed',supporting+tad,G())
        # Discriminating incorrect compositions and map omissions.
        naive=integrate(pull0(nu,nl,u,m2,lam,yy,masses,vc).order(1))+reduced
        reject('apply_original_normalizer_directly_to_reduced_packet',naive,pulled_c)
        wrongbridge=add_shifts(vc,vs,1)
        wrong=integrate(pull0(nu,nl,u,m2,lam,yy,masses,wrongbridge).order(1))+reduced
        reject('use_sum_instead_of_difference_bridge',wrong,pulled_c)
        for label,kw in (('extra_Majorana_half',{'sterile_factor':F(1,2)}),('omit_y_charge_conjugate',{'drop_yc':True}),('commute_retained_mass_through_flavor',{'wrong_mass_order':True})):
            vw=shifts(nu,nl,u,m2,lam,ff,yy,masses,'S',**kw)
            wrong=integrate(pull0(nu,nl,u,m2,lam,yy,masses,vw).order(1))+original
            reject(label,wrong,reduced)
        stationary=(m2+2*lam*u*u==0)
        if stationary: check('stationary_limit','gradient correction vanishes',tad,G())
        else:
            reject('omit_offshell_scalar_tadpole',only_j_s+original,reduced)
            bs_bad=(bridge[0],bridge[1],-bridge[2],-bridge[3])
            reject('wrong_scalar_bridge_sign',integrate(pull0(nu,nl,u,m2,lam,yy,masses,bs_bad).order(1))+reduced,pulled_c)
        trials.append({'p':[str(x) for x in p],'p_squared':str(sum(ETA[mu]*p[mu]**2 for mu in range(4))),
                       'm_H_squared':str(m2),'scalar_gradient_factor':str(m2+2*lam*u*u),'stationary':stationary,
                       'real_fermion_components':n,'original_kernel_monomials':len(original.a),
                       'reduced_kernel_monomials':len(reduced.a),'normalized_kernel_monomials':len(pulled_c.a),
                       'nonlinear_tadpole_monomials':len(tad.a)})
    result={'status':'PASS','arithmetic':'exact Gaussian rational Fourier exterior algebra, standard library only',
            'positive_checks':sum(checks.values()),'checks_by_category':checks,
            'negative_checks':sum(mutations.values()),'mutation_rejections':mutations,
            'total_checks':sum(checks.values())+sum(mutations.values()),'neutral_trials':trials,
            'all_field_bridge':'Xi_S=Xi_C+kappa(v_C-v_S), with maps evaluated on Xi_C',
            'new_direct_substitution_scope':'neutral fermionic quadratic Hessian of declared neutrino/Yukawa/Majorana/scalar action, with nonlinear scalar shift and independent probes',
            'limitations':['full interacting bridge uses prior checked residual and first variation, not a new component expansion of every electroweak vertex',
                           'first order in selected kappa F and inverse-scale order three','no gauge fixing or ghost pullback',
                           'no ordinary quantum 1PI coordinate-covariance claim','no new loops, complete physical basis, finite amplitude, or beta function']}
    Path(__file__).with_name('FIELD_MAP_COMPOSITION_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Exact action-composition replay completed in %.3f seconds' % (time.perf_counter()-started),flush=True)
    print(json.dumps(result,indent=2))


if __name__=='__main__': run()
