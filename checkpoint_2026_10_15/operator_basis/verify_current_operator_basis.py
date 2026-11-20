#!/usr/bin/env python3
"""Four-dimensional current-sector quotient using literal exterior components.

No repository imports or computer-algebra package. Weak and Lorentz invariant
dimensions come from separate infinitesimal-generator linear systems. Actual
identical barred-lepton variables, not distinct placeholders, enter the
Grassmann algebra and its coefficient-vector rank calculation.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json


class Q:
    __slots__=('r','i')
    def __init__(self,r=0,i=0):
        if isinstance(r,Q): self.r,self.i=r.r,r.i
        else: self.r,self.i=F(r),F(i)
    def __add__(self,o):
        o=Q(o); return Q(self.r+o.r,self.i+o.i)
    __radd__=__add__
    def __neg__(self): return Q(-self.r,-self.i)
    def __sub__(self,o): return self+-Q(o)
    def __rsub__(self,o): return Q(o)+-self
    def __mul__(self,o):
        o=Q(o); return Q(self.r*o.r-self.i*o.i,self.r*o.i+self.i*o.r)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=Q(o); d=o.r*o.r+o.i*o.i
        if not d: raise ZeroDivisionError
        return Q((self.r*o.r+self.i*o.i)/d,(self.i*o.r-self.r*o.i)/d)
    def __bool__(self): return bool(self.r or self.i)
    def __eq__(self,o):
        o=Q(o); return self.r==o.r and self.i==o.i


I=Q(0,1); ETA=(1,-1,-1,-1)
ONE=((Q(1),Q()),(Q(),Q(1)))
PAULI=(((Q(),Q(1)),(Q(1),Q())),((Q(),-I),(I,Q())),((Q(1),Q()),(Q(),Q(-1))))
T=tuple(tuple(tuple(x*F(1,2) for x in row) for row in m) for m in PAULI)
SIGMA=(ONE,)+PAULI
BAR=(ONE,)+tuple(tuple(tuple(-x for x in row) for row in m) for m in PAULI)
SLOTS=tuple(product(range(2),repeat=4))
LL={s:sum((ETA[mu]*BAR[mu][s[0]][s[1]]*BAR[mu][s[2]][s[3]] for mu in range(4)),Q()) for s in SLOTS}
LR={s:sum((ETA[mu]*BAR[mu][s[0]][s[1]]*SIGMA[mu][s[2]][s[3]] for mu in range(4)),Q()) for s in SLOTS}
W0={s:Q(s[0]==s[1] and s[2]==s[3]) for s in SLOTS}
WX={s:Q(s[0]==s[3] and s[2]==s[1]) for s in SLOTS}


def clean(v): return {k:Q(x) for k,x in v.items() if x}
def plus(a,b,scale=1):
    out=dict(a)
    for k,v in b.items(): out[k]=out.get(k,Q())+scale*v
    return clean(out)
def times(c,a): return clean({k:c*v for k,v in a.items()})
def rank(vectors):
    pivots={}
    for vector in vectors:
        v=clean(vector)
        while v:
            p=min(v)
            if p not in pivots:
                pivots[p]={k:x/v[p] for k,x in v.items()}
                break
            v=plus(v,pivots[p],-v[p])
    return len(pivots)


def ward_rows(actions):
    """Rows of the action on coefficient tensors; actions may be independent."""
    rows=[]
    for group in actions:
        for generator in T:
            for target in SLOTS:
                row={}
                for slot,transpose,sign in group:
                    for x in range(2):
                        original=list(target); original[slot]=x; original=tuple(original)
                        c=generator[x][target[slot]] if transpose else generator[target[slot]][x]
                        row[original]=row.get(original,Q())+sign*c
                rows.append(clean(row))
    return rows


def contract_row(row,tensor): return sum((v*tensor.get(k,Q()) for k,v in row.items()),Q())


def exterior(word):
    """Canonical monomial and its actual anticommuting permutation sign."""
    if len(set(word))!=len(word): return (),0
    parity=sum(word[i]>word[j] for i in range(len(word)) for j in range(i+1,len(word)))
    return tuple(sorted(word)),(-1 if parity%2 else 1)


def term(out,phi,word,c):
    if not c: return
    monomial,sign=exterior(word)
    if not sign: return
    key=(phi,monomial); out[key]=out.get(key,Q())+c*sign


# Species codes make all fermion components independent unless deliberately
# aliased: N=0, L=1, Lbar=2, Q=3, Qbar=4, fR=5, fRbar=6.
def lbar(a,w,s): return (2,a,w,s)
def llep(a,w,s): return (1,a,w,s)
def sterile(i,s): return (0,i,0,s)
def fbar(species,b,w,s): return lbar(b,w,s) if species=='L' else (4 if species=='Q' else 6,b,w,s)
def ffield(species,g,w,s): return llep(g,w,s) if species=='L' else (3 if species=='Q' else 5,g,w,s)


def operator(species,a,b,g,i,exchange=False,scalar=False):
    out={}; right=species not in ('L','Q')
    for wa in range(2):
        for wc in range(1 if right else 2):
            phi=wc if exchange else wa
            wd=wa if exchange else wc
            for sd,sn,sb,sf in SLOTS:
                if scalar:
                    c=Q(sd==sf and sb==sn)
                    word=(lbar(a,wa,sd),ffield(species,g,wd,sf),fbar(species,b,wc,sb),sterile(i,sn))
                else:
                    c=(LR if right else LL)[sd,sn,sb,sf]
                    word=(lbar(a,wa,sd),sterile(i,sn),fbar(species,b,wc,sb),ffield(species,g,wd,sf))
                term(out,phi,word,c)
    return clean(out)


def source_direct(species,nl,nf,nn,ff,g2,b2,yf,normalize=1,yl=-F(1,2)):
    """Literal j JY and sum_I jI JI before weak completeness or Fierz."""
    out={}; right=species not in ('L','Q')
    for a,b,i in product(range(nl),range(nf),range(nn)):
        for wa,wp in product(range(2),repeat=2):
            for wc in range(1 if right else 2):
                for wd in range(1 if right else 2):
                    weak=Q(-b2*yl*yf*F(1,6)*(wa==wp)*(wc==wd))
                    if not right:
                        weak+=sum((-g2*F(1,6)*normalize*t[wa][wp]*t[wc][wd] for t in T),Q())
                    for sd,sn,sb,sf in SLOTS:
                        c=weak*ff[a,i]*(LR if right else LL)[sd,sn,sb,sf]
                        word=(lbar(a,wa,sd),sterile(i,sn),fbar(species,b,wc,sb),ffield(species,b,wd,sf))
                        term(out,wp,word,c)
    return clean(out)


def mapped(species,nl,nf,nn,ff,g2,b2,yf,flip_exchange=False,symmetrize=False):
    out={}; right=species not in ('L','Q')
    for a,b,g,i in product(range(nl),range(nf),range(nf),range(nn)):
        if species=='L':
            c=(g2-b2)*F(1,24)*ff[a,i]*(b==g)-g2*F(1,12)*ff[b,i]*(a==g)*(-1 if flip_exchange else 1)
            if symmetrize:
                c=((g2-b2)*F(1,24)*(ff[a,i]*(b==g)+ff[b,i]*(a==g))
                   -g2*F(1,12)*(ff[b,i]*(a==g)+ff[a,i]*(b==g)))*F(1,2)
            out=plus(out,operator(species,a,b,g,i),c)
        elif right:
            out=plus(out,operator(species,a,b,g,i),b2*yf*F(1,12)*ff[a,i]*(b==g))
        else:
            c0=(g2*F(1,24)+b2*yf*F(1,12))*ff[a,i]*(b==g)
            cx=-g2*F(1,12)*ff[a,i]*(b==g)
            out=plus(out,operator(species,a,b,g,i),c0)
            out=plus(out,operator(species,a,b,g,i,exchange=True),cx)
    return out


checks={}; rejects={}; ranks=[]
def check(cat,label,a,b):
    if a!=b: raise RuntimeError('FAILED '+cat+': '+label)
    checks[cat]=checks.get(cat,0)+1
def reject(label,a,b):
    if a==b: raise RuntimeError('MUTATION SURVIVED '+label)
    rejects[label]=rejects.get(label,0)+1


def run():
    # The two weak invariants and the unique chiral Lorentz tensors are
    # determined by separate full tensor-space invariant equations.
    ww=ward_rows((((0,False,1),(1,True,-1),(2,False,1),(3,True,-1)),))
    lw=ward_rows((((0,False,1),(2,False,1)),((1,True,1),(3,True,1))))
    rw=ward_rows((((0,False,1),(3,True,-1)),((1,True,1),(2,False,-1))))
    check('invariant_dimensions','weak kernel on all sixteen components',16-rank(ww),2)
    check('invariant_dimensions','LL Lorentz kernel on all sixteen components',16-rank(lw),1)
    check('invariant_dimensions','LR Lorentz kernel on all sixteen components',16-rank(rw),1)
    check('invariant_dimensions','weak candidate independence',rank((W0,WX)),2)
    for name,rows,tensors in (('weak',ww,(W0,WX)),('LL',lw,(LL,)),('LR',rw,(LR,))):
        for tensor in tensors:
            for row in rows: check('full_component_invariance',name,contract_row(row,tensor),Q())
    for s in SLOTS:
        a,b,c,d=s
        actual=sum((t[a][b]*t[c][d] for t in T),Q())
        check('weak_completeness','ordered generator tensor',actual,WX[s]*F(1,2)-W0[s]*F(1,4))
    # Chiral Fierz, including all anticommuting reordering signs.
    for nl,nn in ((1,1),(2,2),(3,1)):
        oo=[]; xx=[]
        for a,b,g,i in product(range(nl),range(nl),range(nl),range(nn)):
            op=operator('L',a,b,g,i); ex=operator('L',a,b,g,i,exchange=True)
            check('identical_LL_exterior_fierz','all weak/spinor components',ex,operator('L',b,a,g,i))
            reject('wrong_identical_LL_Fierz_sign',ex,times(-1,operator('L',b,a,g,i)))
            oo.append(op); xx.append(ex)
        expected=nl**3*nn
        check('identical_LL_rank','ordered basis rank',rank(oo),expected)
        check('identical_LL_rank','no extra exchange directions',rank(oo+xx),expected)
        ranks.append({'species':'L','n_L':nl,'n_N':nn,'columns':len(oo),'rank':expected,'rank_with_exchange':expected})
    # Distinct Q fields have two independent weak contractions per flavor slot.
    for nl,nf,nn in ((1,1,1),(2,2,2)):
        basis=[operator('Q',a,b,g,i,exchange=x) for a,b,g,i,x in product(range(nl),range(nf),range(nf),range(nn),(False,True))]
        expected=2*nl*nf**2*nn
        check('distinct_Q_rank','all weak/spinor components',rank(basis),expected)
        ranks.append({'species':'Q','n_L':nl,'n_f':nf,'n_N':nn,'columns':len(basis),'rank':expected})
    for species in ('u','d','e'):
        basis=[]
        for a,b,g,i in product(range(2),range(2),range(2),range(2)):
            op=operator(species,a,b,g,i); sc=operator(species,a,b,g,i,scalar=True)
            check('right_current_exterior_fierz',species,op,times(-2,sc))
            reject('wrong_right_current_scalar_sign',op,times(2,sc))
            basis.append(op)
        check('right_current_rank',species,rank(basis),16)
        ranks.append({'species':species,'n_L':2,'n_f':2,'n_N':2,'columns':16,'rank':16})
    # Direct current substitutions, tested separately in the independent g^2
    # and g'^2 coefficient directions, then with mixed rational couplings.
    for species,yf,nl,nf,nn in (('L',-F(1,2),3,3,2),('Q',F(1,6),2,2,2),
                               ('u',F(2,3),2,2,2),('d',-F(1,3),2,2,2),('e',-1,2,2,2)):
        ff={(a,i):Q((a+1)*(i+2)-1,2*a-i+1) for a,i in product(range(nl),range(nn))}
        for g2,b2 in ((1,0),(0,1),(F(7,3),F(5,2))):
            direct=source_direct(species,nl,nf,nn,ff,g2,b2,yf)
            expected=mapped(species,nl,nf,nn,ff,g2,b2,yf)
            check('representative_coefficient_map',species,direct,expected)
            if b2: reject('wrong_hypercharge_sign_'+species,direct,source_direct(species,nl,nf,nn,ff,g2,b2,yf,yl=F(1,2)))
            if g2 and species in ('L','Q'):
                reject('Pauli_instead_of_T_'+species,direct,source_direct(species,nl,nf,nn,ff,g2,b2,yf,normalize=4))
            if g2 and species=='L':
                reject('wrong_L_exchange_coefficient_sign',direct,mapped(species,nl,nf,nn,ff,g2,b2,yf,flip_exchange=True))
                reject('spurious_L_barred_flavor_symmetry',direct,mapped(species,nl,nf,nn,ff,g2,b2,yf,symmetrize=True))
    ff={(0,0):Q(2,3)}
    lone=source_direct('L',1,1,1,ff,F(7,3),F(5,2),-F(1,2))
    check('one_flavor_L_limit','nonzero ordered representative',lone,times(-(F(7,3)+F(5,2))*F(1,24)*ff[0,0],operator('L',0,0,0,0)))
    reject('one_flavor_L_incorrectly_zero',lone,{})
    # An explicit obstruction to imposing symmetry of the two barred flavors.
    reject('O0_incorrectly_symmetric_in_barred_flavors',operator('L',0,1,0,0),operator('L',1,0,0,0))
    result={'status':'PASS','arithmetic':'exact Gaussian rational; standard library only',
            'scope':'algebraic four-dimensional nonderivative four-fermion+phi sectors of the selected current representative',
            'positive_checks':sum(checks.values()),'positive_checks_by_category':checks,
            'negative_checks':sum(rejects.values()),'mutation_rejections':rejects,
            'total_checks':sum(checks.values())+sum(rejects.values()),'component_ranks':ranks,
            'invariant_dimensions':{'weak_doublet':2,'LL_Lorentz':1,'LR_Lorentz':1},
            'L_basis_dimension':'n_L^3 n_N','Q_basis_dimension':'2 n_L n_Q^2 n_N',
            'singlet_basis_dimension_each':'n_L n_f^2 n_N',
            'L_exchange_identity':'Oex(alpha,beta,gamma,i)=O0(beta,alpha,gamma,i)',
            'L_coefficient':'[(g^2-gprime^2) F_alpha,i delta_beta,gamma /24 - g^2 F_beta,i delta_alpha,gamma /12]',
            'one_flavor_L_coefficient':'-(g^2+gprime^2)F/24',
            'Q_coefficients':{'O0':'(3g^2+gprime^2) F delta /72','Oex':'-g^2 F delta /12'},
            'right_current_vector_coefficients':{'u':'gprime^2 F delta /18','d':'-gprime^2 F delta /36','e':'-gprime^2 F delta /12'},
            'right_scalar_fierz_factor':-2,
            'limitations':['selected source representative only','no derivative/EOM or other-operator quotient',
                           'strictly four-dimensional Lorentz and Fierz algebra','no finite evanescent prescription',
                           'no new loops, beta function, observable, or novelty claim']}
    Path(__file__).with_name('CURRENT_OPERATOR_BASIS_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': run()
