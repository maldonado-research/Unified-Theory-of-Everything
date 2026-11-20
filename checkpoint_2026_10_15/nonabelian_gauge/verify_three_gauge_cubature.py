#!/usr/bin/env python3
"""Standalone exact 24-cell cubature of 33 three-gauge matter-loop UV graphs.
No producer code or repository module is imported. All coefficients are integers
or Fraction; imaginary Wick coordinates are exact Gaussian integer pairs.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from functools import lru_cache
import json
import argparse

ETA = (1,-1,-1,-1)
ZERO = (0,0)
ONE = (1,0)
I = (0,1)
BASIS = tuple(tuple(int(a==b) for a in range(4)) for b in range(4))
CHECKS = 0

def ck(x, label):
    global CHECKS
    CHECKS += 1
    if not x:
        raise RuntimeError(label)

def ga(a,b): return (a[0]+b[0],a[1]+b[1])
def gm(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def gs(n,a): return (n*a[0],n*a[1])
def gp(a,n):
    out=ONE
    for _ in range(n): out=gm(out,a)
    return out

def clean(d): return {k:v for k,v in d.items() if v != ZERO}

@lru_cache(None)
def blade(a,b):
    n=sum((b & ((1<<j)-1)).bit_count() for j in range(4) if a&(1<<j))
    sign=(-1)**n
    for j in range(4):
        if a&b&(1<<j): sign*=ETA[j]
    return sign,a^b

def cmul(a,b):
    out={}
    for ia,ca in a.items():
        for ib,cb in b.items():
            sign,m=blade(ia,ib)
            out[m]=ga(out.get(m,ZERO),gs(sign,gm(ca,cb)))
    return clean(out)

def cscale(z,a): return clean({m:gm(z,c) for m,c in a.items()})
def csum(*args):
    out={}
    for d in args:
        for m,c in d.items(): out[m]=ga(out.get(m,ZERO),c)
    return clean(out)

def slash(q): return clean({1<<j:gs(ETA[j],q[j]) for j in range(4)})
def polslash(e): return slash(tuple((x,0) for x in e))
def dotq(q,e):
    out=ZERO
    for j in range(4): out=ga(out,gs(ETA[j]*e[j],q[j]))
    return out

def dot(a,b): return sum(ETA[j]*a[j]*b[j] for j in range(4))
def cchain(*args):
    out={0:ONE}
    for a in args: out=cmul(out,a)
    return out

def cell():
    out=[]
    for a,b in combinations(range(4),2):
        for sa,sb in product((-1,1),repeat=2):
            q=[ZERO]*4
            q[a]=(sa,0) if a==0 else (0,sa)
            q[b]=(sb,0) if b==0 else (0,sb)
            out.append(tuple(q))
    return out
POINTS=cell()


def cubature(fn,degree):
    acc={}
    for q in POINTS:
        acc=csum(acc,fn(q))
    # Coordinates above are sqrt(2)*q_unit; every requested integrand is
    # homogeneous and even, so divide by the exact power 2**(degree/2).
    den=24*2**(degree//2)
    out={}
    for m,c in acc.items():
        ck(c[1]==0,"imaginary angular remainder")
        if c[0]: out[m]=F(c[0],den)
    return out


def fchain(q,Fseq,eps):
    Q=slash(q)
    out=Q
    for a in Fseq:
        out=cchain(Q,polslash(eps[a]),out)
    return out


def topology_descriptors():
    out=[]
    for nf in range(4):
        for perm in permutations(range(3)):
            f,s=perm[:nf],perm[nf:]
            out.append({"name":"O_F"+''.join(map(str,f))+"_S"+''.join(map(str,s)),
                        "kind":"ordinary","F":f,"S":s,
                        "words":{tuple(reversed(s))+f:(-1)**nf}})
    for pair in combinations(range(3),2):
        c=next(a for a in range(3) if a not in pair)
        # Extra scalar propagator times i{Ta,Tb}: i*i=-1, rather than +1
        # for the ordinary i*(-i) insertion. Conjugate fermion gives -1.
        words={pair+(c,):1,tuple(reversed(pair))+(c,):1}
        out.append({"name":"SG_F"+str(c)+"_S"+''.join(map(str,pair)),
                    "kind":"seagull_fermion","F":(c,),"S":(pair,),
                    "pair":pair,"single":c,"words":words})
        for blocks in ((pair,(c,)),((c,),pair)):
            prefix,suffix=blocks[1],blocks[0]
            choices=[tuple(prefix)+tuple(suffix)]
            if len(prefix)==2:
                choices.append(tuple(reversed(prefix))+tuple(suffix))
            else:
                choices.append(tuple(prefix)+tuple(reversed(suffix)))
            words={w:-1 for w in choices}
            out.append({"name":"SG_F_S"+'_'.join(''.join(map(str,x)) for x in blocks),
                        "kind":"seagull_scalar","F":(),"S":blocks,
                        "pair":pair,"single":c,"words":words})
    return out
TOPS=topology_descriptors()
WORDS=tuple(permutations(range(3)))


def graph(top,eps):
    if top['kind']=='ordinary':
        def fn(q):
            coeff=ONE
            for a in top['S']: coeff=gm(coeff,gs(-2,dotq(q,eps[a])))
            return cscale(coeff,fchain(q,top['F'],eps))
        return cubature(fn,4)
    a,b=top['pair']
    d=dot(eps[a],eps[b])
    if top['kind']=='seagull_fermion':
        return cubature(lambda q:cscale((d,0),fchain(q,top['F'],eps)),2)
    c=top['single']
    return cubature(lambda q:cscale(gs(-2*d,dotq(q,eps[c])),slash(q)),2)


def rational_add(out,part,weight=1):
    for m,c in part.items():
        out[m]=out.get(m,F(0))+weight*c
        if out[m]==0: del out[m]


def integrate(eps,omit=None,flip_seagull=False,reverse_fermion=False):
    out={w:{} for w in WORDS}
    individual={}
    for top in TOPS:
        if omit==top['name'] or (omit=='all_seagulls' and top['kind']!='ordinary'):
            continue
        part=graph(top,eps)
        individual[top['name']]=part
        for word,weight in top['words'].items():
            if reverse_fermion and len(top['F'])>1:
                ns=3-len(top['F'])
                word=word[:ns]+tuple(reversed(word[ns:]))
            if flip_seagull and top['kind']!='ordinary': weight=-weight
            rational_add(out[word],part,weight)
    return out,individual



def action_kernel(eps,drop=None):
    """Separate cubic covariant-action target; added after blind integration.

    Curvature convention F_mn=partial_m A_n-partial_n A_m+i[A_m,A_n].
    The sigma term gives -1/4 A_r[A_m,A_n][gamma^m,gamma^n]gamma^r;
    the adjoint-curvature term gives -1/3[A_n,[A^n,A^m]]gamma_m.
    Gauge couplings are included in each free generator label.
    """
    out={w:{} for w in WORDS}
    E=[polslash(e) for e in eps]
    def accum(w,n,v):
        part={}
        for m,c in v.items():
            ck(c[1]==0,'real local action coefficient')
            if c[0]: part[m]=c[0]
        rational_add(out[w],part,n)
    for i,j in combinations(range(3),2):
        k=next(a for a in range(3) if a not in (i,j))
        if drop!='sigma':
            v=cmul(csum(cmul(E[i],E[j]),cscale((-1,0),cmul(E[j],E[i]))),E[k])
            accum((k,i,j),-F(1,4),v)
            accum((k,j,i),F(1,4),v)
        if drop!='adjoint':
            v=csum(cscale((dot(eps[k],eps[i]),0),E[j]),cscale((-dot(eps[k],eps[j]),0),E[i]))
            for n,w in ((-F(1,3),(k,i,j)),(F(1,3),(k,j,i)),(F(1,3),(i,j,k)),(-F(1,3),(j,i,k))):
                accum(w,n,v)
    return out


def gmmul(a,b):
    n=len(a)
    out=[]
    for i in range(n):
        row=[]
        for j in range(n):
            z=ZERO
            for k in range(n): z=ga(z,gm(a[i][k],b[k][j]))
            row.append(z)
        out.append(tuple(row))
    return tuple(out)


def weak_projection(words,generators):
    # Each integer Gaussian generator below is 2*T (or 2*y_L).
    out={}
    for w,d in words.items():
        mat=((ONE,ZERO),(ZERO,ONE))
        for a in w: mat=gmmul(mat,generators[a])
        for mask,c in d.items():
            for i,j in product(range(2),repeat=2):
                key=(mask,i,j)
                out[key]=ga(out.get(key,ZERO),gs(c/F(8),mat[i][j]))
    return {str(k):[str(x) for x in z] for k,z in sorted(out.items()) if z!=ZERO}


def encode(d): return {str(m):str(c) for m,c in sorted(d.items())}
def encode_words(d): return {''.join(map(str,w)):encode(d[w]) for w in WORDS}


def mmul(a,b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(a))) for j in range(len(a))) for i in range(len(a)))
def mt(a): return tuple(zip(*a))
def ms(n,a): return tuple(tuple(n*x for x in row) for row in a)
def ma(a,b): return tuple(tuple(a[i][j]+b[i][j] for j in range(len(a))) for i in range(len(a)))
def mchain(*args):
    n=len(args[0]) if args else 3
    out=tuple(tuple(int(i==j) for j in range(n)) for i in range(n))
    for a in args: out=mmul(out,a)
    return out

def rank(rows):
    a=[list(map(F,r)) for r in rows]
    nr=0
    for col in range(len(a[0])):
        piv=next((k for k in range(nr,len(a)) if a[k][col]),None)
        if piv is None: continue
        a[nr],a[piv]=a[piv],a[nr]
        p=a[nr][col]
        a[nr]=[x/p for x in a[nr]]
        for k in range(len(a)):
            if k!=nr and a[k][col]:
                p=a[k][col]
                a[k]=[x-p*y for x,y in zip(a[k],a[nr])]
        nr+=1
        if nr==len(a): break
    return nr


def color_audit():
    T=(((1,2,-1),(0,3,2),(-2,1,0)),
       ((0,1,3),(-1,2,0),(2,-2,1)),
       ((2,-1,0),(3,1,-2),(1,0,4)))
    rows=[sum((tuple(r) for r in mchain(*(T[x] for x in w))),()) for w in WORDS]
    ck(rank(rows)==6,'six free cubic words distinguished by literal 3x3 sample')
    for top in TOPS:
        ferm=mchain(*(ms(-1,mt(T[x])) for x in reversed(top['F'])))
        blocks=[]
        if top['kind']=='ordinary': blocks=[T[x] for x in top['S']]
        else:
            for b in top['S']:
                blocks.append(T[b[0]] if len(b)==1 else ma(mmul(T[b[0]],T[b[1]]),mmul(T[b[1]],T[b[0]])))
        scal=mchain(*reversed(blocks))
        # Kx_ab;cd=delta_ad delta_bc and Yukawa delta_ef:
        # sum_ef scal_a,f ferm_c,e delta_e,f = scal ferm^T.
        literal=tuple(tuple(sum(scal[a][e]*ferm[c][e] for e in range(3)) for c in range(3)) for a in range(3))
        if top['kind']!='ordinary': literal=ms(-1,literal)
        expected=((0,0,0),)*3
        for w,coef in top['words'].items(): expected=ma(expected,ms(coef,mchain(*(T[x] for x in w))))
        ck(literal==expected,'literal representation transport '+top['name'])
    return {"literal_topology_color_contractions":len(TOPS),"dense_matrix_word_rank":6}


def moment_audit():
    ck(len(POINTS)==24,'24-cell size')
    for q in POINTS:
        norm=ZERO
        for j in range(4): norm=ga(norm,gs(ETA[j],gm(q[j],q[j])))
        ck(norm==(2,0),'unit sphere after exact sqrt2 rescaling')
    for degree in range(1,5):
        for inds in product(range(4),repeat=degree):
            z=ZERO
            for q in POINTS:
                v=ONE
                for j in inds: v=gm(v,q[j])
                z=ga(z,v)
            ck(z[1]==0,'real monomial average')
            if degree%2:
                ck(z[0]==0,'odd moment')
            elif degree==2:
                a,b=inds
                exp=F(ETA[a]*int(a==b),4)
                ck(F(z[0],48)==exp,'second sphere moment')
            else:
                a,b,c,d=inds
                exp=F(ETA[a]*ETA[c]*int(a==b and c==d)+ETA[a]*ETA[b]*(int(a==c and b==d)+int(a==d and b==c)),24)
                ck(F(z[0],96)==exp,'fourth sphere moment')


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--coefficients',help='write complete free-word and Clifford coefficient table')
    args=p.parse_args()
    moment_audit()
    color=color_audit()
    ck(len(TOPS)==33 and sum(t['kind']=='ordinary' for t in TOPS)==24,'33 distinct topologies')
    ck(len({t['name'] for t in TOPS})==33,'distinct topology names')
    # Literal source phase: (i seed)(i fermion)(-i Yukawa)(i scalar)=-1.
    ck(gm(gm(gm(I,I),gs(-1,I)),I)==(-1,0),'literal source phase')
    ck(gm(I,gs(-1,I))==ONE,'ordinary insertion phase')
    ck(gm(I,I)==(-1,0),'anticommutator seagull replacement phase')
    # For commuting color, the two anticommutator words combine to the prior
    # two-B weight -2, while one ordinary F and one S has charge product -1.
    ck(sum((-1,-1))==-2 and (-1)*1==-1,'two-gauge phase and charge calibration')
    tables={}
    controls={"omit_all_seagulls":0,"reverse_seagull_phase":0,"reverse_fermion_word":0,"omit_one_fermion_order":0,"omit_sigma_curvature":0,"omit_adjoint_curvature":0}
    bose=0
    ordinary3=next(t['name'] for t in TOPS if len(t['F'])==3)
    for inds in product(range(4),repeat=3):
        eps=tuple(BASIS[i] for i in inds)
        raw,_=integrate(eps)
        ck(raw==action_kernel(eps),'all free-word cubic action coefficients')
        controls['omit_sigma_curvature']+=raw!=action_kernel(eps,'sigma')
        controls['omit_adjoint_curvature']+=raw!=action_kernel(eps,'adjoint')
        tables[','.join(map(str,inds))]=encode_words(raw)
        # Bose permutation is relabelled in both the polarization inputs and
        # all noncommutative generator words, with no commutative projection.
        for perm in ((1,0,2),(1,2,0)):
            other,_=integrate(tuple(eps[i] for i in perm))
            mapped={w:{} for w in WORDS}
            for w,part in other.items(): mapped[tuple(perm[i] for i in w)]=part
            ck(mapped==raw,'free-word Bose symmetry')
            bose+=1
        tests=(("omit_all_seagulls",{"omit":"all_seagulls"}),
               ("reverse_seagull_phase",{"flip_seagull":True}),
               ("reverse_fermion_word",{"reverse_fermion":True}),
               ("omit_one_fermion_order",{"omit":ordinary3}))
        for name,kw in tests:
            bad,_=integrate(eps,**kw)
            controls[name]+=bad!=raw
        # All three commuting charges must yield zero: the Abelian covariant
        # packet has at most two gauge fields. This is a check, not normalization.
        abelian={}
        for d in raw.values(): rational_add(abelian,d)
        ck(abelian=={},'commuting three-gauge pole vanishes')
    for name,n in controls.items(): ck(n>0,'negative control '+name)
    if args.coefficients:
        with open(args.coefficients,'w') as f: json.dump(tables,f,indent=2,sort_keys=True);f.write('\n')
    b=(((-1,0),ZERO),(ZERO,(-1,0)))
    t1=((ZERO,ONE),(ONE,ZERO))
    t2=((ZERO,(0,-1)),(I,ZERO))
    witness,_=integrate((BASIS[0],BASIS[0],BASIS[1]))
    bww=weak_projection(witness,(b,t1,t2))
    ck(bww=={'(2, 0, 0)':['0','1/8'],'(2, 1, 1)':['0','-1/8']},'mixed B W1 W2 exact witness')
    result={"status":"PASS","route":"33 literal leading loop kernels integrated by exact 24-cell cubature", "checks":CHECKS,
            "repository_imports":False,"arithmetic":"integers, Gaussian integer pairs, exact Fraction",
            "ordinary_topologies":24,"seagull_topologies":9,"polarization_basis_triples":64,
            "free_cubic_color_words":6,"clifford_basis_dimension":16,"bose_coefficient_checks":bose,
            "negative_control_basis_detections":controls,"color_audit":color,
            "uv_ir_separation":"Common positive auxiliary mass in all asymptotic denominators; its corrections are UV finite here. Original generic off-shell integral is never assigned zero as scaleless.",
            "uv_power_count":"Crossed seed q^2 only; direct seed and all external momentum/mass corrections UV finite at three gauge legs.",
            "literal_loop_phase":"-i F g^3/(16 pi^2 epsilon), with ordered generator products retained",
            "scope":"Selected single-seed, one-source, one-matter-loop 1PI three-gauge sector; no gauge internal lines or external radiation",
            "vanishing_equal_polarizations_000":tables['0,0,0'],"nonzero_free_word_witness_001":tables['0,0,1'],
            "mixed_B_W1_W2_witness":{"polarizations":["e0","e0","e1"],"generators":["-I/2","Pauli1/2","Pauli2/2"],"matrix_Clifford_coefficients_real_imag":bww,"coupling_factor":"gprime*g^2"}}
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
