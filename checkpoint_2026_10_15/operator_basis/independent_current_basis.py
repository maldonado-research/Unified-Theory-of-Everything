#!/usr/bin/env python3
"""Independent exact component/exterior-algebra audit of selected current descendants.

Python standard library only. No repository code or stored coefficients are
imported. Four-dimensional Weyl matrices and anticommuting field components
are expanded before identities or flavor ranks are compared.
"""
from fractions import Fraction as R
from functools import lru_cache
from itertools import product
from pathlib import Path
import json


def gc(z): return z if isinstance(z, tuple) else (R(z), R(0))
def add(z,w):
    z,w=gc(z),gc(w); return z[0]+w[0],z[1]+w[1]
def mul(z,w):
    z,w=gc(z),gc(w); return z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0]
def scale(c,z): return mul(c,z)
def real(z):
    z=gc(z)
    if z[1]: raise RuntimeError('unexpected imaginary component coefficient')
    return z[0]


I=(R(0),R(1)); ZERO=(R(0),R(0))
eye=((1,0),(0,1))
pauli=(((0,1),(1,0)),((0,scale(-1,I)),(I,0)),((1,0),(0,-1)))
sigma=(eye,)+pauli
barsigma=(eye,)+tuple(tuple(tuple(scale(-1,z) for z in row) for row in mat) for mat in pauli)
eta=(1,-1,-1,-1)
eps=((0,1),(-1,0))
T=tuple(tuple(tuple(scale(R(1,2),z) for z in row) for row in mat) for mat in pauli)


def spin_tensor(right):
    second=sigma if right else barsigma
    out={}
    for a,b,c,d in product(range(2),repeat=4):
        z=ZERO
        for mu in range(4): z=add(z,scale(eta[mu],mul(barsigma[mu][a][b],second[mu][c][d])))
        out[a,b,c,d]=real(z)
    return out


LL=spin_tensor(False); LR=spin_tensor(True)
WT={}
for a,b,c,d in product(range(2),repeat=4):
    z=ZERO
    for j in range(3): z=add(z,mul(T[j][a][b],T[j][c][d]))
    WT[a,b,c,d]=real(z)


def gen(tag,flavor,weak,spin,color=0): return tag,flavor,weak,color,spin
def canonical(word):
    if len(set(word))!=len(word): return None,0
    inversions=sum(word[i]>word[j] for i in range(len(word)) for j in range(i+1,len(word)))
    return tuple(sorted(word)), -1 if inversions%2 else 1
def term(out,phi,word,c):
    key,sign=canonical(word)
    if not sign or not c: return
    key=phi,key
    out[key]=out.get(key,R(0))+sign*c
    if not out[key]: del out[key]
def polyadd(a,b,c=R(1)):
    out=dict(a)
    for k,v in b.items():
        out[k]=out.get(k,R(0))+c*v
        if not out[k]: del out[k]
    return out
def polyscale(c,a): return {k:c*v for k,v in a.items() if c*v}


SPEC={'L':(False,1,R(-1,2)), 'Q':(False,3,R(1,6)),
      'uR':(True,3,R(2,3)), 'dR':(True,3,R(-1,3)), 'eR':(True,1,R(-1))}


@lru_cache(None)
def operator(species,kind,alpha,beta,gamma,i):
    """D direct, E exchanged weak contraction, or W literal Pauli insertion."""
    right,colors,_=SPEC[species]
    if right and kind!='D': raise ValueError('weak singlet has no exchange operator')
    out={}; tensor=LR if right else LL
    if right:
        weak=[(a,a,0,0,R(1)) for a in range(2)]
    elif kind=='D':
        weak=[(a,a,c,c,R(1)) for a,c in product(range(2),repeat=2)]
    elif kind=='E':
        weak=[(a,b,b,a,R(1)) for a,b in product(range(2),repeat=2)]
    elif kind=='W':
        weak=[(a,b,c,d,w) for (a,b,c,d),w in WT.items() if w]
    else: raise ValueError('operator kind')
    bar='Lb' if species=='L' else species+'b'
    field='L' if species=='L' else species
    for a,b,c,d,w in weak:
        for (s,t,u,v),z in tensor.items():
            if not z: continue
            for color in range(colors):
                word=(gen('Lb',alpha,a,s),gen('N',i,0,t),
                      gen(bar,beta,c,u,color),gen(field,gamma,d,v,color))
                term(out,b,word,w*z)
    return out


counts={}; mutations={}
def check(category,label,a,b):
    if a!=b: raise RuntimeError('FAILED '+category+': '+label)
    counts[category]=counts.get(category,0)+1
def reject(label,a,b):
    if a==b: raise RuntimeError('MUTATION SURVIVED '+label)
    mutations[label]=mutations.get(label,0)+1


def pack_add(out,polynomial,gauge,fa,fi,c):
    """Keep each F_(fa fi) and each squared gauge coupling as a formal symbol."""
    for k,v in polynomial.items():
        key=gauge,fa,fi,k
        out[key]=out.get(key,R(0))+c*v
        if not out[key]: del out[key]


def raw_source(species,nl,ns,nf):
    out={}; right,_,yf=SPEC[species]
    for alpha,i,beta in product(range(nl),range(ns),range(nf)):
        d=operator(species,'D',alpha,beta,beta,i)
        pack_add(out,d,'gp2',alpha,i,yf/12)
        if not right: pack_add(out,operator(species,'W',alpha,beta,beta,i),'g2',alpha,i,-R(1,6))
    return out


def reduced_source(species,nl,ns,nf):
    out={}; right,_,yf=SPEC[species]
    for alpha,beta,gamma,i in product(range(nl),range(nf),range(nf),range(ns)):
        d=operator(species,'D',alpha,beta,gamma,i)
        if beta==gamma:
            pack_add(out,d,'gp2',alpha,i,yf/12)
            if not right: pack_add(out,d,'g2',alpha,i,R(1,24))
        if species=='L':
            if alpha==gamma: pack_add(out,d,'g2',beta,i,-R(1,12))
        elif not right and beta==gamma:
            pack_add(out,operator(species,'E',alpha,beta,gamma,i),'g2',alpha,i,-R(1,12))
    return out


def witness(species,kind,alpha,beta,gamma,i):
    if SPEC[species][0]:
        word=(gen('Lb',alpha,0,0),gen('N',i,0,0),gen(species+'b',beta,0,0),gen(species,gamma,0,0))
    else:
        first,second=(0,1) if kind=='D' else (1,0)
        bar='Lb' if species=='L' else species+'b'
        field='L' if species=='L' else species
        word=(gen('Lb',alpha,first,0),gen('N',i,0,0),gen(bar,beta,second,1),gen(field,gamma,1,1))
    canonical_word,sign=canonical(word)
    if not sign: raise RuntimeError('zero witness')
    return (0,canonical_word),sign*2


def rank_certificate(species,nl,ns,nf):
    kinds=('D','E') if species=='Q' else ('D',)
    columns=[(kind,a,b,c,i) for kind,a,b,c,i in product(kinds,range(nl),range(nf),range(nf),range(ns))]
    rows=[witness(species,*col) for col in columns]
    if len(set(k for k,_ in rows))!=len(rows): raise RuntimeError('duplicate witness rows')
    determinant=1
    # Explicit square component matrix is diagonal. This is a rank certificate,
    # with no numerical random points or tolerance and no matrix inversion.
    for column,col in enumerate(columns):
        p=operator(species,*col)
        for row,(key,expected) in enumerate(rows):
            check('explicit_diagonal_rank_matrix',species,p.get(key,R(0)),R(expected) if row==column else R(0))
        determinant*=rows[column][1]
    return {'species':species,'n_L':nl,'n_N':ns,'spectator_flavors':nf,
            'columns':len(columns),'rank':len(columns),'component_minor_determinant':str(determinant),
            'witness_rows':[{'column':list(col),'phi':key[0],
                             'grassmann_generators':[list(g) for g in key[1]],'diagonal':diag}
                            for col,(key,diag) in zip(columns,rows)]}


def run():
    for a,b,c,d in product(range(2),repeat=4):
        check('four_dimensional_spin_components','LL',LL[a,b,c,d],R(2*eps[a][c]*eps[b][d]))
        check('four_dimensional_spin_components','LR',LR[a,b,c,d],R(2*(a==d)*(b==c)))
        check('pauli_components','T=tau/2',WT[a,b,c,d],R((a==d)*(b==c),2)-R((a==b)*(c==d),4))
    nl,ns,nf=3,2,3
    for a,b,c,i in product(range(nl),range(nl),range(nl),range(ns)):
        e=operator('L','E',a,b,c,i); flipped=operator('L','D',b,a,c,i)
        check('grassmann_LL_Fierz','ordered barred-flavor exchange',e,flipped)
        reject('wrong_Fierz_minus',e,polyscale(-1,flipped))
    for species in ('L','Q'):
        for a,b,c,i in product(range(nl),range(nf),range(nf),range(ns)):
            raw=operator(species,'W',a,b,c,i)
            target=polyadd(polyscale(R(1,2),operator(species,'E',a,b,c,i)),operator(species,'D',a,b,c,i),-R(1,4))
            check('full_component_weak_reconstruction',species,raw,target)
    for species in SPEC:
        raw=raw_source(species,nl,ns,nf); reduced=reduced_source(species,nl,ns,nf)
        check('symbolic_F_and_gauge_source_reconstruction',species,raw,reduced)
        reject('omit_'+species+'_current',raw,{})
    # A one-flavor LL channel is nonzero and has both negative squared couplings.
    one=raw_source('L',1,1,1); d=operator('L','D',0,0,0,0); target={}
    pack_add(target,d,'g2',0,0,-R(1,24)); pack_add(target,d,'gp2',0,0,-R(1,24))
    check('one_flavor_limit','negative g2+gp2 coefficient',one,target)
    reject('one_flavor_claimed_zero',d,{})
    reject('extra_identical_lepton_half',one,{k:v/2 for k,v in target.items()})
    qsource=raw_source('Q',nl,ns,nf)
    reject('extra_color_multiplicity',qsource,{k:3*v for k,v in qsource.items()})
    reject('Pauli_instead_of_half_Pauli_generators',qsource,
           {k:(4*v if k[0]=='g2' else v) for k,v in qsource.items()})
    q_direct={}
    for a,i,b in product(range(nl),range(ns),range(nf)):
        dop=operator('Q','D',a,b,b,i)
        pack_add(q_direct,dop,'g2',a,i,R(1,24)); pack_add(q_direct,dop,'gp2',a,i,R(1,72))
    reject('omit_Q_weak_exchange',qsource,q_direct)
    reject('reverse_hypercharge_current_sign',qsource,
           {k:(-v if k[0]=='gp2' else v) for k,v in qsource.items()})
    # Flipping the F index in the exchange descendant changes the actual source.
    wrong={}
    for a,b,c,i in product(range(nl),range(nl),range(nl),range(ns)):
        d=operator('L','D',a,b,c,i)
        if b==c:
            pack_add(wrong,d,'g2',a,i,R(1,24)); pack_add(wrong,d,'gp2',a,i,-R(1,24))
        if a==c: pack_add(wrong,d,'g2',a,i,-R(1,12))
    reject('exchange_keeps_wrong_F_flavor',raw_source('L',nl,ns,nf),wrong)
    certificates=[rank_certificate(species,nl,ns,nf) for species in SPEC]
    result={'status':'PASS','method':'exact rational Weyl/Pauli components, canonical Grassmann monomials, explicit diagonal component minors',
            'tests_by_category':counts,'positive_checks':sum(counts.values()),
            'mutation_rejections':mutations,'negative_checks':sum(mutations.values()),
            'total_checks':sum(counts.values())+sum(mutations.values()),
            'rank_certificates':certificates,
            'general_sector_dimensions':{'L':'n_L^3 n_N','Q':'2 n_L n_N n_Q^2',
                                         'uR':'n_L n_N n_u^2','dR':'n_L n_N n_d^2','eR':'n_L n_N n_e^2'},
            'scope':['four-dimensional nonderivative phi+four-fermion current channel',
                     'only selected gauge-current descendant; coefficients carry common original kappa and h.c.',
                     'ranks refer to unrestricted ordered-flavor component polynomials, not generic nonzero seed coefficients',
                     'not complete EFT EOM closure, dimensional evanescent prescription, loop matching or physical running']}
    target=Path(__file__).with_name('INDEPENDENT_CURRENT_BASIS_RESULTS.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rank_certificates'},indent=2))
    print('Certified component ranks: '+', '.join(x['species']+'='+str(x['rank']) for x in certificates))


if __name__=='__main__': run()
