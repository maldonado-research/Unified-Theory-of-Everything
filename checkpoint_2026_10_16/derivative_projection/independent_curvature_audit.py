#!/usr/bin/env python3
"""Independent exact invariant/component and covariant-jet curvature audit.

Standard library only. No project imports, stored coefficients or removable
assertions. Ranks concern the fixed derivative-on-scalar placement only.
"""
from fractions import Fraction as R
from itertools import combinations, permutations, product
from pathlib import Path
import json


def z(v): return v if isinstance(v,tuple) else (R(v),R(0))
def za(a,b):
    a,b=z(a),z(b); return a[0]+b[0],a[1]+b[1]
def zm(a,b):
    a,b=z(a),z(b); return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def zs(a): return a[0],-a[1]
def nz(a): return bool(a[0] or a[1])
I=(R(0),R(1)); Z=z(0); ETA=(1,-1,-1,-1)
def zsum(values):
    out=Z
    for v in values: out=za(out,v)
    return out
def zero(n,m=None): return [[Z for _ in range(n if m is None else m)] for _ in range(n)]
def eye(n): return [[z(i==j) for j in range(n)] for i in range(n)]
def ma(a,b,c=1): return [[za(v,zm(c,w)) for v,w in zip(row,other)] for row,other in zip(a,b)]
def ms(c,a): return [[zm(c,v) for v in row] for row in a]
def mm(a,b): return [[zsum(zm(a[i][k],b[k][j]) for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]
def comm(a,b): return ma(mm(a,b),mm(b,a),-1)
def dagger_matrix(a): return [[zs(a[j][i]) for j in range(len(a))] for i in range(len(a[0]))]
def mv(a,v): return [zsum(zm(a[i][j],v[j]) for j in range(len(v))) for i in range(len(a))]


def parity(word): return (-1)**sum(word[j]>word[k] for j in range(len(word)) for k in range(j+1,len(word)))
def epsilon(a,b,c,d): return 0 if len({a,b,c,d})!=4 else parity((a,b,c,d))
PAIRS=list(combinations(range(4),2))
P1=[[z(0),z(1)],[z(1),z(0)]]; P2=[[Z,zm(-1,I)],[I,Z]]; P3=[[z(1),Z],[Z,z(-1)]]
PAULI=(P1,P2,P3); T=tuple(ms(R(1,2),p) for p in PAULI)
SIGMA=(eye(2),)+PAULI
BS=(eye(2),)+tuple(ms(-1,p) for p in PAULI)
GAMMA=[]
for mu in range(4):
    a=zero(4)
    for s,t in product(range(2),repeat=2):
        a[s][t+2]=SIGMA[mu][s][t]; a[s+2][t]=BS[mu][s][t]
    GAMMA.append(a)
G5=ms(I,mm(mm(mm(GAMMA[0],GAMMA[1]),GAMMA[2]),GAMMA[3]))
PL=[[z(i==j and i<2) for j in range(4)] for i in range(4)]
SPIN=tuple(tuple(ms(zm(I,R(1,2)),comm(GAMMA[mu],GAMMA[nu])) for nu in range(4)) for mu in range(4))
TRIPLE={}
for mu,nu,rho in product(range(4),repeat=3):
    a=zero(4)
    for order in permutations(range(3)):
        labels=(mu,nu,rho)
        a=ma(a,mm(mm(GAMMA[labels[order[0]]],GAMMA[labels[order[1]]]),GAMMA[labels[order[2]]]),R(parity(order),6))
    TRIPLE[mu,nu,rho]=a


def ec(c): return {():z(c)} if nz(z(c)) else {}
def field(*key): return {(tuple(key),):z(1)}
def ea(a,b,c=1):
    out=dict(a)
    for key,v in b.items():
        out[key]=za(out.get(key,Z),zm(c,v))
        if not nz(out[key]): del out[key]
    return out
def es(c,a): return {key:zm(c,v) for key,v in a.items() if nz(zm(c,v))}
def em(a,b):
    out={}
    for ka,va in a.items():
        for kb,vb in b.items():
            word=ka+kb
            if len(set(word))!=len(word): continue
            key=tuple(sorted(word)); coeff=zm(parity(word),zm(va,vb))
            out[key]=za(out.get(key,Z),coeff)
            if not nz(out[key]): del out[key]
    return out
def etotal(values):
    out={}
    for v in values: out=ea(out,v)
    return out
def row_times(row,a): return [etotal(es(a[s][t],row[s]) for s in range(len(row))) for t in range(len(a[0]))]
def matrix_column(a,col): return [etotal(es(a[s][t],col[t]) for t in range(len(col))) for s in range(len(a))]
def bilinear(row,a,col): return etotal(es(a[s][t],em(row[s],col[t])) for s,t in product(range(4),repeat=2) if nz(a[s][t]))
def contract(row,weak,spin,col,phi):
    return etotal(es(zm(weak[a][b],phi[b]),bilinear(row[a],spin,col))
                  for a,b in product(range(2),repeat=2) if nz(weak[a][b]) and nz(phi[b]))


counts={}; mutations={}
def check(category,label,a,b):
    if a!=b: raise RuntimeError('FAILED '+category+': '+label)
    counts[category]=counts.get(category,0)+1
def reject(label,a,b):
    if a==b: raise RuntimeError('MUTATION SURVIVED '+label)
    mutations[label]=mutations.get(label,0)+1


def sparse_rank(rows):
    pivots={}
    for initial in rows:
        row={k:R(v) for k,v in initial.items() if v}
        while row:
            first=min(row)
            if first not in pivots:
                divisor=row[first]; pivots[first]={k:v/divisor for k,v in row.items()}; break
            coeff=row[first]
            for k,v in pivots[first].items():
                row[k]=row.get(k,R(0))-coeff*v
                if not row[k]: del row[k]
    return len(pivots)


def lorentz_kernel():
    keys=list(product(range(4),range(4),PAIRS)); index={k:j for j,k in enumerate(keys)}
    constraints=[]
    for a,b in PAIRS:
        m=[[0]*4 for _ in range(4)]; m[a][b]=ETA[b]; m[b][a]=-ETA[a]
        rows=[{} for _ in keys]
        for column,(r,s,(u,v)) in enumerate(keys):
            terms={}
            def put(key,c): terms[key]=terms.get(key,0)+c
            for t in range(4):
                put((t,s,(u,v)),m[r][t])             # b contravariant
                put((r,t,(u,v)),-m[t][s])           # Dphi covariant
                for first,second,c in ((t,v,-m[t][u]),(u,t,-m[t][v])):
                    if first==second: continue
                    put((r,s,tuple(sorted((first,second)))),c*(1 if first<second else -1))
            for key,c in terms.items():
                if c: rows[index[key]][column]=c
        constraints.extend(rows)
    rk=sparse_rank(constraints)
    # Explicit ordinary and dual invariant coefficients on the same 96 monomials.
    ordinary=[]; dual=[]
    for r,s,(a,b) in keys:
        ordinary.append(ETA[s]*((r==a and s==b)-(r==b and s==a)))
        dual.append(ETA[r]*epsilon(r,s,a,b))
    for label,vector in (('ordinary',ordinary),('dual',dual)):
        for row in constraints:
            check('Lorentz_kernel_invariant_vectors',label,sum(c*vector[k] for k,c in row.items()),0)
    return {'unknowns':len(keys),'constraint_rows':len(constraints),'rank':rk,'nullity':len(keys)-rk}


def weak_kernels():
    generators=(((0,1),(0,0)),((0,0),(1,0)),((1,0),(0,-1)))
    wb=generators
    result=[]
    for adjoint in (False,True):
        keys=list(product(range(2),range(2),range(3))) if adjoint else list(product(range(2),repeat=2))
        index={k:j for j,k in enumerate(keys)}; constraints=[]
        for gen in generators:
            ad=[[0]*3 for _ in range(3)]
            for j,w in enumerate(wb):
                bracket=[[sum(gen[a][k]*w[k][b]-w[a][k]*gen[k][b] for k in range(2)) for b in range(2)] for a in range(2)]
                ad[0][j],ad[1][j],ad[2][j]=bracket[0][1],bracket[1][0],bracket[0][0]
            rows=[{} for _ in keys]
            for column,key in enumerate(keys):
                a,b=key[:2]; terms={}
                def put(k,c): terms[k]=terms.get(k,0)+c
                for c in range(2):
                    put((c,b)+key[2:],-gen[c][a]); put((a,c)+key[2:],gen[b][c])
                if adjoint:
                    for j in range(3): put((a,b,j),ad[key[2]][j])
                for k,c in terms.items():
                    if c: rows[index[k]][column]=c
            constraints.extend(rows)
        rk=sparse_rank(constraints)
        vector=[wb[key[2]][key[0]][key[1]] if adjoint else int(key[0]==key[1]) for key in keys]
        for row in constraints:
            check('weak_kernel_invariant_vectors','W' if adjoint else 'B',sum(c*vector[k] for k,c in row.items()),0)
        result.append({'strength':'W' if adjoint else 'B','unknowns':len(keys),'constraint_rank':rk,'nullity':len(keys)-rk})
    return result


def component_operator(channel,dual):
    out={}
    for a,b,s,t,nu,(r,q) in product(range(2),range(2),range(2),range(2),range(4),PAIRS):
        for label in range(3) if channel=='W' else range(1):
            weak=T[label][a][b] if channel=='W' else z(a==b)
            coeff=Z
            for mu in range(4):
                c=ETA[mu]*epsilon(mu,nu,r,q) if dual else ETA[nu]*((mu==r and nu==q)-(mu==q and nu==r))
                coeff=za(coeff,zm(c,BS[mu][s][t]))
            value=zm(weak,coeff)
            if nz(value): out[(channel,label,r,q,nu,b,a,s,t)]=value
    return out


def hermitian(seed):
    out=ms(R(seed+1,7),eye(2))
    for j in range(3): out=ma(out,PAULI[j],R(((seed+2)*(j+1))%7-3,j+2))
    return out
def jets(seed):
    aa=[hermitian(seed+mu) for mu in range(4)]
    da=[[hermitian(seed+3*mu+nu+5) for nu in range(4)] for mu in range(4)]
    dda=[[[hermitian(seed+4*(rho+mu)+nu+9) for nu in range(4)] for mu in range(4)] for rho in range(4)]
    ff=[[ma(ma(da[mu][nu],da[nu][mu],-1),comm(aa[mu],aa[nu]),I) for nu in range(4)] for mu in range(4)]
    df=[]
    for rho in range(4):
        layer=[]
        for mu in range(4):
            row=[]
            for nu in range(4):
                value=ma(dda[rho][mu][nu],dda[rho][nu][mu],-1)
                value=ma(value,ma(comm(da[rho][mu],aa[nu]),comm(aa[mu],da[rho][nu])),I)
                row.append(ma(value,comm(aa[rho],ff[mu][nu]),I))
            layer.append(row)
        df.append(layer)
    return ff,df
def dual_lower(ff):
    return [[ms(R(ETA[mu]*ETA[nu],2),
                [[zsum(zm(epsilon(mu,nu,r,s),ff[r][s][a][b]) for r,s in product(range(4),repeat=2)) for b in range(2)] for a in range(2)])
             for nu in range(4)] for mu in range(4)]


def run():
    eps=((0,1),(-1,0))
    expected_c=[[z(eps[s][t] if s<2 and t<2 else -eps[s-2][t-2] if s>=2 and t>=2 else 0)
                 for t in range(4)] for s in range(4)]
    check('charge_conjugation_components','C=i gamma2 gamma0',ms(I,mm(GAMMA[2],GAMMA[0])),expected_c)
    for mu,nu in product(range(4),repeat=2):
        check('Clifford_components','anticommutator',ma(mm(GAMMA[mu],GAMMA[nu]),mm(GAMMA[nu],GAMMA[mu])),ms(2*ETA[mu] if mu==nu else 0,eye(4)))
    for mu,nu,rho in product(range(4),repeat=3):
        rhs=zero(4)
        for s in range(4): rhs=ma(rhs,mm(GAMMA[s],G5),zm(I,epsilon(mu,nu,rho,s)*ETA[s]))
        check('triple_gamma_and_chirality','G=i epsilon gamma_lower gamma5',TRIPLE[mu,nu,rho],rhs)
        check('triple_gamma_and_chirality','gamma5 PL=-PL',mm(G5,PL),ms(-1,PL))
    component_matrix=[[BS[mu][s][t] for mu in range(4)] for s,t in product(range(2),repeat=2)]
    determinant=Z
    for order in permutations(range(4)):
        v=z(parity(order))
        for r in range(4): v=zm(v,component_matrix[r][order[r]])
        determinant=za(determinant,v)
    reject('Weyl_bilinear_not_a_full_vector_module',determinant,Z)
    lorentz=lorentz_kernel(); weak=weak_kernels()
    check('invariant_dimensions','Lorentz two',lorentz['nullity'],2)
    for item in weak: check('invariant_dimensions','weak '+item['strength'],item['nullity'],1)
    operators=[component_operator(c,d) for c in ('B','W') for d in (False,True)]
    rows=[('B',0,0,1,1,0,0,0,0),('B',0,2,3,1,0,0,0,0),
          ('W',2,0,1,1,0,0,0,0),('W',2,2,3,1,0,0,0,0)]
    diagonal=(-1,1,R(-1,2),R(1,2))
    for r,key in enumerate(rows):
        for c,op in enumerate(operators): check('actual_chiral_component_minor',str((r,c)),op.get(key,Z),z(diagonal[r]) if r==c else Z)
    image=[]
    for c in range(2):
        op={}
        for key in set(operators[2*c])|set(operators[2*c+1]):
            value=zm(R(1,2),za(operators[2*c].get(key,Z),zm(I,operators[2*c+1].get(key,Z))))
            if nz(value): op[key]=value
        image.append(op)
    reject('ordinary_and_dual_are_identical',operators[0],operators[1])
    reject('selected_chiral_image_zero',image[0],{})
    for r,key in enumerate((rows[0],rows[2])):
        for c,op in enumerate(image):
            check('selected_image_component_minor',str((r,c)),op.get(key,Z),z((R(-1,2),R(-1,4))[r]) if r==c else Z)
    gp_image=R(3,5); g_image=R(7,6); y_image=R(-1,2)
    raw_f_column=[z(gp_image*y_image/2),zm(I,gp_image*y_image/2),z(g_image/2),zm(I,g_image/2)]
    check('fixed_coupling_raw_F_column','B ordinary dual tie',raw_f_column[1],zm(I,raw_f_column[0]))
    check('fixed_coupling_raw_F_column','W ordinary dual tie',raw_f_column[3],zm(I,raw_f_column[2]))
    check('fixed_coupling_raw_F_column','B W common F tie',zm(g_image,raw_f_column[0]),zm(gp_image*y_image,raw_f_column[2]))
    reject('selected_fixed_coupling_raw_F_column_zero',raw_f_column,[Z]*4)
    for trial in range(2):
        ff,df=jets(trial+1); star=dual_lower(ff)
        check('Hodge_star','star squared minus one',dual_lower(star),[[ms(-1,f) for f in row] for row in ff])
        selected=[[ms(R(1,2),ma(ff[mu][nu],star[mu][nu],I)) for nu in range(4)] for mu in range(4)]
        check('Hodge_star','selected field has eigenvalue minus i',dual_lower(selected),[[ms(zm(-1,I),v) for v in row] for row in selected])
        for rho,mu,nu in product(range(4),repeat=3):
            check('nonAbelian_Bianchi','full adjoint derivative',ma(ma(df[rho][mu][nu],df[mu][nu][rho]),df[nu][rho][mu]),zero(2))
        row=[[field('lb',a,s) if s>=2 else {} for s in range(4)] for a in range(2)]
        drow=[[[field('dlb',rho,a,s) if s>=2 else {} for s in range(4)] for a in range(2)] for rho in range(4)]
        n=[field('n',s) if s<2 else {} for s in range(4)]
        dn=[[field('dn',rho,s) if s<2 else {} for s in range(4)] for rho in range(4)]
        phi=[z((R(2+trial),R(1))),z((R(1),R(-2-trial)))]
        dp=[[z((R(mu+1),R(trial+1))),z((R(2-mu),R(-mu-1)))] for mu in range(4)]
        kr=[[etotal(es(zm(-1,I),row_times(drow[rho][a],GAMMA[rho])[s]) for rho in range(4)) for s in range(4)] for a in range(2)]
        ancol=[etotal(es(I,matrix_column(GAMMA[rho],dn[rho])[s]) for rho in range(4)) for s in range(4)]
        qv={}; qe={}; hstar={}; qr={}; kterm={}; qdiv={}; qan={}; odrow={}; divergence={}
        qvp={}; qep={}; qer={}; bv={}; be={}
        for mu,nu in product(range(4),repeat=2):
            qv=ea(qv,contract(row,ff[mu][nu],GAMMA[mu],n,dp[nu]),-ETA[nu])
            hstar=ea(hstar,contract(row,star[mu][nu],GAMMA[mu],n,dp[nu]),ETA[nu])
            qr=ea(qr,contract(drow[nu],ff[mu][nu],GAMMA[mu],n,phi),-ETA[nu])
            kterm=ea(kterm,contract(kr,ff[mu][nu],SPIN[mu][nu],n,phi))
            qdiv=ea(qdiv,contract(row,df[nu][nu][mu],GAMMA[mu],n,phi),ETA[nu])
            qan=ea(qan,contract(row,ff[mu][nu],SPIN[mu][nu],ancol,phi))
            qvp=ea(qvp,contract(row,ff[mu][nu],GAMMA[mu],dn[nu],phi),ETA[nu])
            bv=ea(bv,etotal((contract(drow[nu],ff[mu][nu],GAMMA[mu],n,phi),
                          contract(row,df[nu][mu][nu],GAMMA[mu],n,phi),
                          contract(row,ff[mu][nu],GAMMA[mu],n,dp[nu]),
                          contract(row,ff[mu][nu],GAMMA[mu],dn[nu],phi))),ETA[nu])
            for rho in range(4):
                sg=mm(SPIN[mu][nu],GAMMA[rho])
                gg=TRIPLE[mu,nu,rho]
                qe=ea(qe,contract(row,ff[mu][nu],TRIPLE[mu,nu,rho],n,dp[rho]),R(-1,2))
                qep=ea(qep,contract(row,ff[mu][nu],gg,dn[rho],phi),R(1,2))
                qer=ea(qer,contract(drow[rho],ff[mu][nu],gg,n,phi),R(-1,2))
                be=ea(be,etotal((contract(drow[rho],ff[mu][nu],gg,n,phi),
                                contract(row,df[rho][mu][nu],gg,n,phi),
                                contract(row,ff[mu][nu],gg,n,dp[rho]),
                                contract(row,ff[mu][nu],gg,dn[rho],phi))),R(1,2))
                first=contract(drow[rho],ff[mu][nu],sg,n,phi)
                odrow=ea(odrow,first,zm(I,R(1,4)))
                divergence=ea(divergence,etotal((first,contract(row,df[rho][mu][nu],sg,n,phi),
                                    contract(row,ff[mu][nu],sg,n,dp[rho]),contract(row,ff[mu][nu],sg,dn[rho],phi))))
        check('dual_curvature_component_identity','QE=-i Hstar',qe,es(zm(-1,I),hstar))
        hplus=es(R(-1,2),ea(qv,qe))
        check('selected_chiral_image','Hplus=(H+iHstar)/2',hplus,es(R(1,2),ea(es(-1,qv),es(I,hstar))))
        od=ea(odrow,es(R(1,3),qdiv))
        check('row_Clifford_route','no IBP',od,etotal((qr,es(R(-1,4),kterm),es(R(1,3),qdiv))))
        check('separate_derivative_placements','row triple gamma',qer,ea(qr,es(R(-1,2),kterm)))
        check('separate_derivative_placements','sterile Clifford',ea(qvp,qep),es(R(-1,2),qan))
        check('separate_derivative_placements','ordinary full product rule',qvp,etotal((qr,qv,qdiv,bv)))
        check('separate_derivative_placements','triple full product rule and Bianchi',qep,etotal((qer,qe,be)))
        check('separate_derivative_placements','combined boundary sign',es(R(1,2),ea(bv,be)),es(zm(-1,zm(I,R(1,4))),divergence))
        expected=etotal((hplus,es(R(-1,6),qdiv),es(R(-1,4),qan),es(zm(I,R(1,4)),divergence)))
        check('covariant_IBP_Bianchi_route','full off-shell boundary term',od,expected)
        alternative=etotal((qr,es(R(-1,4),kterm),es(R(1,2),qdiv),es(R(1,4),qan),es(zm(-1,zm(I,R(1,4))),divergence)))
        check('alternative_derivative_placement','all row sterile and curvature derivatives',hplus,alternative)
        reject('reverse_dual_sign',qe,es(I,hstar))
        reject('drop_boundary_term',hplus,etotal((qr,es(R(-1,4),kterm),es(R(1,2),qdiv),es(R(1,4),qan))))
        reject('omit_gauge_divergence',hplus,ea(alternative,qdiv,R(-1,2)))
        reject('omit_sterile_derivative_term',hplus,ea(alternative,qan,R(-1,4)))
        # Independent nonzero variational source and residual slots.
        sb=[[field('sb',a,s) if s<2 else {} for s in range(4)] for a in range(2)]
        eb=[[field('eb',a,s) if s<2 else {} for s in range(4)] for a in range(2)]
        sr=[field('sr',s) if s>=2 else {} for s in range(4)]
        er=[field('er',s) if s>=2 else {} for s in range(4)]
        nr=[field('nr',s) if s>=2 else {} for s in range(4)]
        ksource=[[ea(sb[a][s],eb[a][s]) for s in range(4)] for a in range(2)]
        rsource=[etotal((es(3,nr[s]),sr[s],er[s])) for s in range(4)]
        left={}; right={}; explicit={}; residual={}
        for mu,nu in product(range(4),repeat=2):
            left=ea(left,contract(ksource,ff[mu][nu],SPIN[mu][nu],n,phi),R(-1,4))
            left=ea(left,contract(row,ff[mu][nu],SPIN[mu][nu],rsource,phi),R(1,4))
            explicit=ea(explicit,contract(sb,ff[mu][nu],SPIN[mu][nu],n,phi),R(-1,4))
            explicit=ea(explicit,contract(row,ff[mu][nu],SPIN[mu][nu],
                            [ea(es(3,nr[s]),sr[s]) for s in range(4)],phi),R(1,4))
            residual=ea(residual,contract(eb,ff[mu][nu],SPIN[mu][nu],n,phi),R(-1,4))
            residual=ea(residual,contract(row,ff[mu][nu],SPIN[mu][nu],er,phi),R(1,4))
        check('nonzero_matter_sources','lepton and sterile residual orientation',left,ea(explicit,residual))
        reject('omit_actual_lepton_sterile_sources',left,residual)
        # H^dagger phi=0, while a weak-strength insertion prevents that zero.
        eps=((0,1),(-1,0)); h=[zsum(zm(-eps[a][b],zs(phi[b])) for b in range(2)) for a in range(2)]
        inner=zsum(zm(zs(h[a]),phi[a]) for a in range(2))
        check('charged_lepton_source','Hdag phi zero',inner,Z)
        inserted=zsum(zm(zs(h[a]),zm(T[0][a][b],phi[b])) for a,b in product(range(2),repeat=2))
        reject('weak_inserted_charged_lepton_source_zero',inserted,Z)
        # Explicit physical Yukawa sources, with the literal C convention.
        bar_nr=[etotal(es(eps[s][t],n[s]) for s in range(2)) if t<2 else {} for t in range(4)]
        nr_phys=[{} if s<2 else etotal(es(-eps[s-2][t],field('nb',t)) for t in range(2)) for s in range(4)]
        erbar=[field('charged_bar',s) if s<2 else {} for s in range(4)]
        yn=z((R(1),R(2))); ye=z((R(2),R(-1)))
        sb_nu=[[es(zm(zs(yn),zs(phi[a])),bar_nr[s]) for s in range(4)] for a in range(2)]
        sb_e=[[es(zm(zs(ye),zs(h[a])),erbar[s]) for s in range(4)] for a in range(2)]
        sb_total=[[ea(sb_nu[a][s],sb_e[a][s]) for s in range(4)] for a in range(2)]
        xbar=[etotal(es(phi[a],row[a][s+2]) for a in range(2)) for s in range(2)]
        xc=[{} if s<2 else etotal(es(-eps[s-2][t],xbar[t]) for t in range(2)) for s in range(4)]
        sr_phys=[es(yn,v) for v in xc]
        actual_source={}; source_parts={}; charged_contact={}; charged_B={}; sterile_contact={}; neutrino_column={}
        for mu,nu in product(range(4),repeat=2):
            f=ff[mu][nu]; spin=SPIN[mu][nu]
            actual_source=ea(actual_source,contract(sb_total,f,spin,n,phi),R(-1,4))
            actual_source=ea(actual_source,contract(row,f,spin,[ea(es(3,nr_phys[s]),sr_phys[s]) for s in range(4)],phi),R(1,4))
            sterile_contact=ea(sterile_contact,contract(sb_nu,f,spin,n,phi),R(-1,4))
            charged_contact=ea(charged_contact,contract(sb_e,f,spin,n,phi),R(-1,4))
            trace_half=zm(R(1,2),za(f[0][0],f[1][1])); pure_b=ms(trace_half,eye(2))
            charged_B=ea(charged_B,contract(sb_e,pure_b,spin,n,phi),R(-1,4))
            source_parts=ea(source_parts,contract(row,f,spin,nr_phys,phi),R(3,4))
            source_parts=ea(source_parts,contract(row,f,spin,xc,phi),zm(R(1,4),yn))
            neutrino_column=ea(neutrino_column,contract(row,f,spin,xc,phi),zm(R(1,4),yn))
        source_parts=etotal((source_parts,sterile_contact,charged_contact))
        check('actual_Yukawa_sources','Majorana mass neutrino and charged source slots',actual_source,source_parts)
        check('actual_Yukawa_sources','one-flavor sterile tensor is zero',sterile_contact,{})
        check('actual_Yukawa_sources','pure B charged source zero',charged_B,{})
        reject('delete_weak_charged_Yukawa_source',charged_contact,{})
        reject('omit_neutrino_sterile_Yukawa_source',actual_source,ea(actual_source,neutrino_column,-1))
        # Gauge divergence source split with actual scalar and quark currents.
        qrow=[[field('qb',a,s) if s>=2 else {} for s in range(4)] for a in range(2)]
        qcol=[[field('q',a,s) if s<2 else {} for s in range(4)] for a in range(2)]
        gp=R(3,5); g=R(7,6); yl=R(-1,2)
        div_source={}; div_explicit={}; div_residual={}
        for mu in range(4):
            scalar_y=zm(zm(I,yl),zsum(za(zm(zs(phi[a]),dp[mu][a]),zm(-1,zm(zs(dp[mu][a]),phi[a]))) for a in range(2)))
            jy=ea(ec(scalar_y),etotal(bilinear(qrow[a],GAMMA[mu],qcol[a]) for a in range(2)),R(ETA[mu],6))
            jw=[]
            for k in range(3):
                scalar_w=zm(I,zsum(zm(T[k][a][b],za(zm(zs(phi[a]),dp[mu][b]),zm(-1,zm(zs(dp[mu][a]),phi[b])))) for a,b in product(range(2),repeat=2)))
                current=etotal(es(T[k][a][b],bilinear(qrow[a],GAMMA[mu],qcol[b])) for a,b in product(range(2),repeat=2))
                jw.append(ea(ec(scalar_w),current,ETA[mu]))
            eb=ec(R(mu+2,7)); ew=[ec(R((mu+1)*(k+2),11)) for k in range(3)]
            for a,b in product(range(2),repeat=2):
                sources=etotal((es(gp*gp*yl*(a==b),jy),etotal(es(zm(g*g,T[k][a][b]),jw[k]) for k in range(3))))
                residuals=etotal((es(gp*yl*(a==b),eb),etotal(es(zm(g,T[k][a][b]),ew[k]) for k in range(3))))
                outer=es(phi[b],bilinear(row[a],GAMMA[mu],n))
                div_source=ea(div_source,em(outer,ea(sources,residuals)))
                div_explicit=ea(div_explicit,em(outer,sources))
                div_residual=ea(div_residual,em(outer,residuals))
        check('gauge_source_certificate','scalar quark currents and nonzero gauge residuals',div_source,ea(div_explicit,div_residual))
        reject('replace_gauge_equation_by_free_equation',div_source,div_residual)
        reject('discard_gauge_residual_off_shell',div_source,div_explicit)
    result={'status':'PASS','method':'exact rational complex Clifford/Weyl components, exterior fermions, non-Abelian connection jets and independently solved invariant kernels',
            'tests_by_category':counts,'positive_checks':sum(counts.values()),
            'mutation_rejections':mutations,'negative_checks':sum(mutations.values()),
            'total_checks':sum(counts.values())+sum(mutations.values()),
            'Lorentz_kernel':lorentz,'weak_kernels':weak,
            'fixed_scalar_derivative_placement_rank':4,'independent_channel_chiral_projector_rank':2,
            'fixed_coupling_raw_F_rank_per_flavor':1,
            'fixed_coupling_raw_F_column_real_imaginary_strings':[[str(c) for c in v] for v in raw_f_column],
            'actual_chiral_component_minor_diagonal':['-1','1','-1/2','1/2'],
            'conventions':{'metric':'+---','epsilon_contravariant_0123':1,'gamma5':'i gamma0 gamma1 gamma2 gamma3',
                           'dual_lower':'eta_mu eta_nu epsilon^(mu nu rho sigma) F_(rho sigma)/2',
                           'selected_curvature':'(H(F)+i H(star F))/2; star eigenvalue -i'},
            'scope':['fixed barL N_L Dphi B/W field-strength placement algebraic invariant rank',
                     'alternative derivative placement retains row, sterile, gauge-divergence and source terms',
                     'not full EFT IBP/EOM independence, finite matching, physical beta function or observable']}
    Path(__file__).with_name('INDEPENDENT_CURVATURE_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': run()
