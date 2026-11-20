#!/usr/bin/env python3
"""Exact independent algebraic tensor-sector quotient; SymPy 1.14.0.

No project imports. Run python -I -B -O independent_tensor_basis.py.
Prints a JSON receipt. This does not compute a derivative-EOM EFT quotient.
"""
import itertools
import json
import sys
import sympy as s

I, R = s.I, s.Rational
pauli = [s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
T = [a/2 for a in pauli]
checks, mutations, ranks = [], [], {}

def iszero(x):
    return all(s.expand(a)==0 for a in x) if isinstance(x,s.MatrixBase) else s.expand(x)==0

def check(name, x):
    if not iszero(x):
        raise RuntimeError(name+': exact nonzero residual '+str(x))
    checks.append(name)

def reject(name,x):
    if iszero(x):
        raise RuntimeError(name+': false statement not rejected')
    mutations.append(name)

# Solve all infinitesimal weak-invariance conditions on the *entire* tensor
# coefficient spaces, rather than testing only proposed operators.
l = s.symbols('l0:2')
r = s.symbols('r0:2')
h = s.symbols('h0:2')
w = s.symbols('w0:3')
variables = l+r+h+w
def lie(poly,k):
    dl = I*T[k].T*s.Matrix(l)
    dr = I*T[k].T*s.Matrix(r)
    dh = -I*T[k]*s.Matrix(h)
    dw = [sum(s.LeviCivita(k,j,a)*w[j] for j in range(3)) for a in range(3)]
    variations = list(dl)+list(dr)+list(dh)+dw
    return s.expand(sum(s.diff(poly,v)*d for v,d in zip(variables,variations)))

def invariant_space(include_w):
    mons = [l[a]*r[b]*h[c]*h[d]*(w[k] if include_w else 1)
            for a in range(2) for b in range(2)
            for c in range(2) for d in range(c,2)
            for k in (range(3) if include_w else [0])]
    def coeff(p):
        q=s.Poly(s.expand(p),*variables)
        return s.Matrix([q.coeff_monomial(m) for m in mons])
    blocks = [s.Matrix.hstack(*(coeff(lie(m,k)) for m in mons)) for k in range(3)]
    constraints = s.Matrix.vstack(*blocks)
    null = constraints.nullspace()
    swap = s.Matrix.hstack(*(coeff(m.xreplace(dict(zip(l+r,r+l)))) for m in mons))
    return mons,coeff,constraints,null,swap

phi=s.Matrix(h)
P0=phi*phi.T
P=[t*P0 for t in T]
B=sum(l[a]*r[b]*P0[a,b] for a in range(2) for b in range(2))
W=sum(l[a]*r[b]*P[k][a,b]*w[k] for a in range(2) for b in range(2) for k in range(3))
Wswap=W.xreplace(dict(zip(l+r,r+l)))
for sector,poly,expected in [('B',B,1),('W',W,2)]:
    mons,co,con,null,swap=invariant_space(sector=='W')
    check(sector+' complete weak-invariant nullity',len(null)-expected)
    check(sector+' raw invariant',con*co(poly))
    if sector=='B':
        check('B weak exchange symmetric',swap*co(poly)-co(poly))
        reject('B weak tensor falsely antisymmetric',swap*co(poly)+co(poly))
        rank=s.Matrix.hstack(co(poly)).rank()
    else:
        sym,anti=co((poly+Wswap)/2),co((poly-Wswap)/2)
        check('W symmetric weak intertwiner',con*sym)
        check('W antisymmetric weak intertwiner',con*anti)
        check('W symmetric exchange eigenvalue',swap*sym-sym)
        check('W antisymmetric exchange eigenvalue',swap*anti+anti)
        rank=s.Matrix.hstack(sym,anti).rank()
        reject('W only one weak intertwiner',rank-1)
    check(sector+' displayed invariants span whole nullspace',rank-expected)
    ranks[sector+'_weak']={'ambient':len(mons),'constraint_rank':con.rank(),'nullity':len(null)}

# Explicit chiral Clifford matrices establish the spin exchange sign and
# that six two-form components carry only three independent chiral tensors.
z=s.zeros(2)
g0=s.BlockMatrix([[z,s.eye(2)],[s.eye(2),z]]).as_explicit()
gam=[g0]+[s.BlockMatrix([[z,p],[-p,z]]).as_explicit() for p in pauli]
g5=I*gam[0]*gam[1]*gam[2]*gam[3]
check('chiral gamma5 convention',g5-s.diag(-1,-1,1,1))
C=I*gam[2]*gam[0]
spin=[]
for mu in range(4):
    for nu in range(mu+1,4):
        block=(I*(gam[mu]*gam[nu]-gam[nu]*gam[mu])*C/2)[2:4,2:4]
        check(f'chiral sigma C symmetric {mu}{nu}',block-block.T)
        spin.append(block)
spinrank=s.Matrix.hstack(*(s.Matrix(list(b)) for b in spin)).rank()
check('chiral field-strength spin rank',spinrank-3)
reject('ordinary and dual tensors falsely doubled',spinrank-6)
ranks['chiral_spin']=spinrank
spinmons=[l[a]*r[b]*w[k] for a in range(2) for b in range(2) for k in range(3)]
def spincoeff(poly):
    p=s.Poly(s.expand(poly),*variables)
    return s.Matrix([p.coeff_monomial(m) for m in spinmons])
spinconstraints=s.Matrix.vstack(*[s.Matrix.hstack(*(spincoeff(lie(m,k)) for m in spinmons)) for k in range(3)])
check('unique chiral spin-one Lorentz intertwiner',12-spinconstraints.rank()-1)
epsilon=s.Matrix([[0,1],[-1,0]])
spinpoly=sum(l[a]*r[b]*(T[k]*epsilon)[a,b]*w[k] for a in range(2) for b in range(2) for k in range(3))
check('displayed chiral spin-one intertwiner invariant',spinconstraints*spincoeff(spinpoly))
check('chiral spin-one intertwiner exchange symmetric',spincoeff(spinpoly-spinpoly.xreplace(dict(zip(l+r,r+l)))))
ranks['spin_intertwiner']={'ambient':12,'constraint_rank':spinconstraints.rank(),'nullity':1}

# Sparse polynomials in genuine exterior generators bar L_(flavor,weak,spin),
# commuting scalar monomials, and independent gauge/tensor components.
# Using all six physical tensor components retains their exact Clifford
# relation. No numerical scalar background is used for the independence rank.
def clean(p):
    return {k:s.expand(v) for k,v in p.items() if s.expand(v)!=0}

def add(*terms):
    out={}
    for factor,p in terms:
        for k,v in p.items():
            out[k]=out.get(k,0)+factor*v
    return clean(out)

def raw(alpha,beta,sector):
    out={}
    mats=[s.eye(2)] if sector=='B' else T
    for gauge,mat in enumerate(mats):
        for a,b,c in itertools.product(range(2),repeat=3):
            if mat[a,c]==0:
                continue
            scal=tuple(sorted((c,b)))
            for q,S in enumerate(spin):
                for u,v in itertools.product(range(2),repeat=2):
                    left,right=(alpha,a,u),(beta,b,v)
                    if left==right or S[u,v]==0:
                        continue
                    sign=1 if left<right else -1
                    key=(tuple(sorted((left,right))),scal,sector,gauge,q)
                    out[key]=out.get(key,0)+sign*mat[a,c]*S[u,v]
    return clean(out)

def poly_check(name,p):
    check(name,s.Matrix(list(p.values())) if p else 0)

def poly_reject(name,p):
    reject(name,s.Matrix(list(p.values())) if p else 0)

def poly_rank(polys):
    if not polys:
        return 0
    keys=sorted(set().union(*(set(p) for p in polys)))
    return s.Matrix([[p.get(k,0) for p in polys] for k in keys]).rank()

for n in range(1,5):
    br,wa,ws=[],[],[]
    for a in range(n):
        for b in range(a,n):
            rb,rb_rev=raw(a,b,'B'),raw(b,a,'B')
            poly_check(f'n={n} B flavor exchange {a}{b}',add((1,rb),(1,rb_rev)))
            rw,rw_rev=raw(a,b,'W'),raw(b,a,'W')
            sy=add((R(1,2),rw),(R(1,2),rw_rev))
            an=add((R(1,2),rw),(-R(1,2),rw_rev))
            ws.append(sy)
            if a<b:
                br.append(rb)
                wa.append(an)
            else:
                poly_reject(f'n={n} diagonal charged W source discarded {a}',sy)
    actual={'B':poly_rank(br),'W_A':poly_rank(wa),'W_S':poly_rank(ws),'combined':poly_rank(br+wa+ws)}
    expected={'B':n*(n-1)//2,'W_A':n*(n-1)//2,'W_S':n*(n+1)//2,'combined':n*n+n*(n-1)//2}
    for key in expected:
        check(f'n={n} complete Grassmann coefficient rank {key}',actual[key]-expected[key])
    ranks[str(n)]=actual
    # Generic nonsymmetric complex flavor data, including nonzero diagonals.
    K=s.Matrix(n,n,lambda a,b: (a+1)*(b+2)+I*(2*a-b+1))
    KA,KS=(K-K.T)/2,(K+K.T)/2
    sourceB=add(*[(K[a,b],raw(a,b,'B')) for a in range(n) for b in range(n)])
    sourceW=add(*[(K[a,b],raw(a,b,'W')) for a in range(n) for b in range(n)])
    targetB=add(*[(2*KA[a,b],raw(a,b,'B')) for a in range(n) for b in range(a+1,n)])
    targetW=add(*[(2*KA[a,b],add((R(1,2),raw(a,b,'W')),(-R(1,2),raw(b,a,'W')))) for a in range(n) for b in range(a+1,n)],
                *[(2*KS[a,b],add((R(1,2),raw(a,b,'W')),(R(1,2),raw(b,a,'W')))) for a in range(n) for b in range(a+1,n)],
                *[(K[a,a],raw(a,a,'W')) for a in range(n)])
    poly_check(f'n={n} B independent-pair normalization',add((1,sourceB),(-1,targetB)))
    poly_check(f'n={n} W independent-pair normalization',add((1,sourceW),(-1,targetW)))
    if n>1:
        poly_reject(f'n={n} lost B pair factor two',add((1,sourceB),(-R(1,2),targetB)))
    poly_reject(f'n={n} lost W diagonal',add(*[(K[a,a],raw(a,a,'W')) for a in range(n)]))
    # Independent K units test source-image rank, rather than presuming K arbitrary.
    image=[add((1,raw(a,b,'B')),(1,raw(a,b,'W'))) for a in range(n) for b in range(n)]
    check(f'n={n} arbitrary-K source image rank',poly_rank(image)-n*n)

# Background projection: retain only phi_0^2 and inspect weak/spin slots.
def background(p):
    return {k:v for k,v in p.items() if k[1]==(0,0)}

def fixed_weak(alpha,beta,a,b,sector,gauge):
    out={}
    for q,S in enumerate(spin):
        for u,v in itertools.product(range(2),repeat=2):
            left,right=(alpha,a,u),(beta,b,v)
            if left==right or S[u,v]==0:
                continue
            key=(tuple(sorted((left,right))),(0,0),sector,gauge,q)
            out[key]=out.get(key,0)+(1 if left<right else -1)*S[u,v]
    return clean(out)

for a,b in [(0,0),(0,1),(1,0)]:
    neutralB=fixed_weak(a,b,0,0,'B',0)
    neutralW=fixed_weak(a,b,0,0,'W',2)
    charged1=fixed_weak(a,b,1,0,'W',0)
    charged2=fixed_weak(a,b,1,0,'W',1)
    predicted=add((R(1,2),neutralW),(R(1,2),charged1),(I/2,charged2))
    poly_check(f'neutral B projection {a}{b}',add((1,background(raw(a,b,'B'))),(-1,neutralB)))
    poly_check(f'neutral and charged W projection {a}{b}',add((1,background(raw(a,b,'W'))),(-1,predicted)))
    poly_reject(f'wrong charged W phase {a}{b}',add((1,background(raw(a,b,'W'))),(-R(1,2),neutralW),(-R(1,2),charged1),(I/2,charged2)))

# Selected-seed matrix order, including a symmetric seed and a Hermitian H.
Cseed=s.Matrix([[1,2+I],[2+I,3]])
Y=s.Matrix([[1,I,2],[2,1-I,-I]])
H=Y*Y.conjugate().T
F=Cseed*Y.conjugate()
K=F*Y.T
check('selected K matrix order',K-Cseed*H.T)
check('selected K_A formula',(K-K.T)/2-(Cseed*H.T-H*Cseed)/2)
check('selected K_S formula',(K+K.T)/2-(Cseed*H.T+H*Cseed)/2)
reject('selected flavor matrices commuted',K-H.T*Cseed)
check('source prefactor B independent pair',-R(1,12)*2+R(1,6))
reject('source prefactor incorrectly uses isolated packet',-R(1,12)+R(1,4))
result={'status':'pass','python':sys.version.split()[0],'sympy':s.__version__,
        'scope':'Complete algebraic nonderivative two-lepton two-scalar single-field-strength tensor monomial sector; no derivative-EOM EFT quotient, loop, beta function, or observable.',
        'exact_checks':len(checks),'rejected_mutations':len(mutations),'ranks':ranks,
        'checks':checks,'mutations':mutations}
print(json.dumps(result,indent=2))
