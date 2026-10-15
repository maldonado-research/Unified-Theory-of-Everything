#!/usr/bin/env python3
"""Selected three-external-gauge 1PI leading-q UV residues.

Scoped calculation, exact rational Clifford algebra and word-valued color.
No original massless graph is assigned a value by scaleless integration.
Only the local radial logarithm is taken after homogeneous large-q expansion,
and a common auxiliary mass separates IR from the UV logarithm.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations, product
from math import prod
import json
import random
import argparse

ETA=(1,-1,-1,-1)
Z=(0,0,0,0)
B=tuple(tuple(int(i==j) for i in range(4)) for j in range(4))

def clean(a): return {k:v for k,v in a.items() if v}
def add(*xs):
    out={}
    for a in xs:
        for k,v in a.items(): out[k]=out.get(k,0)+v
    return clean(out)
def scale(c,a): return clean({k:c*v for k,v in a.items()})
@lru_cache(None)
def blade(a,b):
    inv=sum((b&((1<<i)-1)).bit_count() for i in range(4) if a&(1<<i))
    return (-1)**inv*prod(ETA[i] for i in range(4) if a&b&(1<<i)),a^b
def mul(*xs):
    out={0:F(1)}
    for a in xs:
        nxt={}
        for m,c in out.items():
            for n,d in a.items():
                s,k=blade(m,n);nxt[k]=nxt.get(k,0)+s*c*d
        out=clean(nxt)
    return out
def slash(v): return {1<<i:ETA[i]*v[i] for i in range(4) if v[i]}
def dot(a,b): return sum(ETA[i]*a[i]*b[i] for i in range(4))
def const(a): return {(Z,m):v for m,v in a.items()}
def qdot(v): return {(B[i],0):ETA[i]*v[i] for i in range(4) if v[i]}
QS={(B[i],1<<i):ETA[i] for i in range(4)}
Q2={(tuple(2*x for x in B[i]),0):ETA[i] for i in range(4)}
def pmul(*xs):
    out={(Z,0):F(1)}
    for a in xs:
        nxt={}
        for (e,m),c in out.items():
            for (f,n),d in a.items():
                s,k=blade(m,n); key=(tuple(x+y for x,y in zip(e,f)),k)
                nxt[key]=nxt.get(key,0)+s*c*d
        out=clean(nxt)
    return out
@lru_cache(None)
def moment(e):
    if any(x%2 for x in e): return F(0)
    r=sum(e)//2
    return F(prod(ETA[i]**(e[i]//2)*prod(range(1,e[i],2)) for i in range(4)),
             prod(4+2*j for j in range(r)))
def angular(a):
    out={}
    for (e,m),c in a.items():out[m]=out.get(m,0)+c*moment(e)
    return clean(out)
def fermion_chain(order,es):
    # order is Yukawa -> seed; open chain is seed -> Yukawa.
    xs=[QS]
    for i in reversed(order):xs.extend((const(slash(es[i])),QS))
    return pmul(*xs)
def ordinary(order_f,order_s,es,uncancelled=False):
    xs=[fermion_chain(order_f,es)]
    xs.extend(scale(-2,qdot(es[i])) for i in order_s)
    if uncancelled:xs.append(Q2)
    return angular(pmul(*xs))
def seagull(order_f,order_s,pair,es,uncancelled=False):
    xs=[fermion_chain(order_f,es)]
    xs.extend(scale(-2,qdot(es[i])) for i in order_s)
    if uncancelled:xs.append(Q2)
    return scale(dot(es[pair[0]],es[pair[1]]),angular(pmul(*xs)))

def raw_diagrams(es):
    out=[]
    # Four distributions of 3 single vertices: 6 each, not 3! times binomial.
    for nf in range(4):
        for perm in permutations(range(3)):
            fs,ss=perm[:nf],perm[nf:]
            out.append((('O',fs,ss),ordinary(fs,ss,es)))
    # A scalar seagull with scalar ordinary insertion has two line orderings.
    for pair in ((0,1),(0,2),(1,2)):
        k=next(i for i in range(3) if i not in pair)
        out.append((('GS',pair,k,'GthenS'),seagull((),(k,),pair,es)))
        out.append((('GS',pair,k,'SthenG'),seagull((),(k,),pair,es)))
        out.append((('GF',pair,k),seagull((k,),(),pair,es)))
    return out

def word_mul(*ws): return tuple(i for w in ws for i in w)
def color_terms(label,mode=0):
    """Literal transport convention, checked against the two-weak-field audit.

    Return (weight,word), with physical seagull relative factor -1.
    Scalar order is reversed in seed->Yukawa; fermion transpose reverses it.
    """
    if label[0]=='O':
        _,fs,ss=label
        if mode==0:word=word_mul(tuple(reversed(ss)),fs)
        else:word=word_mul(fs,tuple(reversed(ss)))
        return [((-1)**len(fs),word)]
    pair,k=label[1:3]
    gp=(pair,tuple(reversed(pair)))
    if label[0]=='GF':
        if mode==0:return [(1,word_mul(g,(k,))) for g in gp]
        return [(1,word_mul((k,),g)) for g in gp]
    if label[3]=='GthenS': return [(-1,word_mul((k,),g)) for g in gp]
    return [(-1,word_mul(g,(k,))) for g in gp]

def colored_total(es,mode=0,mutation=None):
    out={}
    for label,v in raw_diagrams(es):
        if mutation=='drop_seagulls' and label[0]!='O':continue
        if mutation=='drop_FFF012' and label==('O',(0,1,2),()):continue
        for c,w in color_terms(label,mode):
            if mutation=='wrong_seagull_sign' and label[0]!='O':c=-c
            if mutation=='same_internal_generators' and label[0]=='O':c*=(-1)**len(label[1])
            for mask,val in v.items():
                key=(w,mask);out[key]=out.get(key,0)+c*val
    return clean(out)

def closed_ordinary(fs,ss,es):
    """Derived by QS E QS=2(q.e)QS-q²E; no target action used."""
    E=[slash(e) for e in es]
    if len(fs)==0:
        a,b,c=ss
        return scale(-F(1,3),add(scale(dot(es[a],es[b]),E[c]),
                    scale(dot(es[a],es[c]),E[b]),scale(dot(es[b],es[c]),E[a])))
    if len(fs)==1:
        a=fs[0];b,c=ss
        return scale(F(1,3),add(scale(dot(es[a],es[b]),E[c]),
                    scale(dot(es[a],es[c]),E[b]),scale(-2*dot(es[b],es[c]),E[a])))
    if len(fs)==2:
        a,b=fs;c=ss[0]
        return add(scale(F(1,2),mul(E[b],E[a],E[c])),
                   scale(-F(1,3),add(scale(dot(es[a],es[b]),E[c]),
                    scale(dot(es[a],es[c]),E[b]),scale(-2*dot(es[b],es[c]),E[a]))))
    a,b,c=fs
    return scale(F(1,3),add(scale(dot(es[a],es[b]),E[c]),
                    scale(-2*dot(es[a],es[c]),E[b]),scale(dot(es[b],es[c]),E[a])))

def closed_diagram(label,es):
    if label[0]=='O':return closed_ordinary(label[1],label[2],es)
    pair,k=label[1:3]
    return scale(-F(1,2)*dot(es[pair[0]],es[pair[1]]),slash(es[k]))

def action_target(es,drop=None):
    out={}
    def accum(c,word,v):
        for m,a in v.items():out[(word,m)]=out.get((word,m),0)+c*a
    for i,j in ((0,1),(0,2),(1,2)):
        k=next(x for x in range(3) if x not in (i,j))
        ec=add(mul(slash(es[i]),slash(es[j])),scale(-1,mul(slash(es[j]),slash(es[i]))))
        v=mul(ec,slash(es[k]))
        if drop!='sigma':
            accum(-F(1,4),(k,i,j),v);accum(F(1,4),(k,j,i),v)
        v=add(scale(dot(es[k],es[i]),slash(es[j])),scale(-dot(es[k],es[j]),slash(es[i])))
        if drop!='adjoint':
            for c,w in ((-F(1,3),(k,i,j)),(F(1,3),(k,j,i)),
                        (F(1,3),(i,j,k)),(-F(1,3),(j,i,k))):accum(c,w,v)
    return clean(out)

def encode(x):return {str(k):str(v) for k,v in sorted(x.items())}

CHECKS={}
def check(test,label):
    CHECKS[label]=CHECKS.get(label,0)+1
    if not test:raise RuntimeError('FAILED: '+label)
def rename(a,permutation):
    return clean({(tuple(permutation[i] for i in word),mask):v for (word,mask),v in a.items()})
def abelian(a):
    out={}
    for (_,m),v in a.items():out[m]=out.get(m,0)+v
    return clean(out)
def neutral_charged_projection(a):
    # e=1,g=sqrt(2), so X_A=diag(0,-1), X_+=T_+, X_-=T_-.
    # Evaluation uses loop-derived formal color coefficients, not action_target.
    matrices=(((0,0),(0,-1)),((0,1),(0,0)),((0,0),(1,0)))
    out={}
    for (word,mask),v in a.items():
        m=((1,0),(0,1))
        for label in word:
            n=matrices[label]
            m=tuple(tuple(sum(m[i][k]*n[k][j] for k in range(2)) for j in range(2)) for i in range(2))
        out[mask]=out.get(mask,0)+v*m[0][0]
    return clean(out)
def expected_neutral_charged(es):
    a,p,m=es;EA,EP,EM=map(slash,es)
    return add(scale(F(1,4),mul(add(mul(EA,EM),scale(-1,mul(EM,EA))),EP)),
               scale(F(1,3),add(scale(dot(a,m),EP),scale(dot(a,p),EM),scale(-2*dot(p,m),EA))))
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--samples',type=int,default=12)
    args=parser.parse_args()
    mutations={name:0 for name in ('drop_seagulls','wrong_seagull_sign','drop_FFF012',
       'same_internal_generators','reverse_scalar_fermion_color','drop_sigma_curvature','drop_adjoint_curvature')}
    for a,b in product(range(4),repeat=2):
        check(add(mul(slash(B[a]),slash(B[b])),mul(slash(B[b]),slash(B[a])))
                    ==clean({0:2*ETA[a]*int(a==b)}),'Clifford basis')
    check((1j)*(1j)*(-1j)*(1j)==-1,'baseline source phase')
    check(1j*(-1j)==1,'ordinary insertion phase')
    check(1j*(1j)==-1,'seagull propagator vertex phase before anticommutator')
    for n in (1,2,3):check(angular(pmul(*(Q2 for _ in range(n))))=={0:F(1)},'angular radial trace')
    rng=random.Random(0x3A20261015)
    samples=list(product(B,repeat=3))
    samples.extend(tuple(tuple(F(rng.randint(-5,5),rng.randint(1,3)) for _ in range(4))
                    for _ in range(3)) for _ in range(args.samples))
    witness=None
    for sample,es in enumerate(samples):
        diagrams=raw_diagrams(es)
        check(len(diagrams)==33 and len({label for label,_ in diagrams})==33,'33 distinct topologies')
        for label,raw in diagrams:
            check(raw==closed_diagram(label,es),'individual angular vs analytic')
            if label[0]=='O':
                uncancelled=ordinary(label[1],label[2],es,True)
            else:
                pair,k=label[1:3]
                uncancelled=seagull((k,) if label[0]=='GF' else (),
                              (k,) if label[0]=='GS' else (),pair,es,True)
            check(raw==uncancelled,'cancelled vs uncancelled seed q²')
        raw=colored_total(es);target=action_target(es)
        check(raw==target,'formal free-color-word cubic action match')
        check(abelian(raw)=={},'pure Abelian cubic residue zero')
        check(neutral_charged_projection(raw)==expected_neutral_charged(es),'direct loop A Wplus Wminus neutral projection')
        for perm in permutations(range(3)):
            check(rename(colored_total(tuple(es[i] for i in perm)),perm)==raw,'full Bose symmetry')
        for name in mutations:
            if name=='reverse_scalar_fermion_color':bad=colored_total(es,1)
            elif name=='drop_sigma_curvature':bad=action_target(es,'sigma')
            elif name=='drop_adjoint_curvature':bad=action_target(es,'adjoint')
            else:bad=colored_total(es,mutation=name)
            mutations[name]+=bad!=target
        if sample==1:witness={'polarizations':es,'loop_sum':encode(raw),'action':encode(target),
                          'individual_diagrams':[(str(l),encode(v)) for l,v in diagrams]}
    for name,count in mutations.items():check(count>0,'mutation rejected: '+name)
    print(json.dumps({'status':'PASS','arithmetic':'exact rational Cl(1,3) and formal free color words',
            'route':'actual homogeneous large-q logarithmic angular integration with auxiliary IR mass',
            'checks':CHECKS,'total_checks':sum(CHECKS.values()),
            'positive_checks':sum(v for k,v in CHECKS.items() if not k.startswith('mutation rejected:')),
            'negative_control_rejections':sum(v for k,v in CHECKS.items() if k.startswith('mutation rejected:')),
            'complete_polarization_basis_points':64,'generic_rational_points':args.samples,
            'matter_topologies_per_point':33,'mutations_detected_at_points':mutations,
            'imports_repository_code':False,'result':'33 selected three-gauge 1PI matter diagrams match independently derived cubic action kernel',
            'full_loop_phase':'-i F/(16*pi^2*epsilon) times generator-coupling-valued residue',
            'scope':'single derivative seed, conjugate Yukawa, three external gauge fields; no extra internal gauge propagators or 1PR tree attachments; no physical RGE or finite amplitude',
            'first_nonzero_exact_witness':witness},indent=2,default=str))
if __name__=='__main__':main()
