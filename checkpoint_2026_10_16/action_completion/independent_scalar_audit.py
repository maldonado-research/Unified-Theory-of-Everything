#!/usr/bin/env python3
"""Independent weak-invariant, scalar-jet and Grassmann-component audit.

Exact rational arithmetic, standard library only, no project imports. This
checks only the scalar-current part of the previously sourced packet.
"""
from fractions import Fraction as R
from itertools import product, combinations_with_replacement
from pathlib import Path
import json


def z(x): return x if isinstance(x,tuple) else (R(x),R(0))
def za(x,y):
    x,y=z(x),z(y); return x[0]+y[0],x[1]+y[1]
def zm(x,y):
    x,y=z(x),z(y); return x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def zs(x): return x[0],-x[1]
def nz(x): return bool(x[0] or x[1])
I=(R(0),R(1)); Z=(R(0),R(0))


def constant(c): return {((),()):z(c)} if nz(z(c)) else {}
def bos(*key): return {((tuple(key),),()):z(1)}
def odd(*key): return {((),(tuple(key),)):z(1)}
def add(a,b,c=1):
    out=dict(a)
    for key,v in b.items():
        out[key]=za(out.get(key,Z),zm(c,v))
        if not nz(out[key]): del out[key]
    return out
def scale(c,a): return {key:zm(c,v) for key,v in a.items() if nz(zm(c,v))}
def total(items):
    out={}
    for a in items: out=add(out,a)
    return out
def mul(a,b):
    out={}
    for (ca,ga),va in a.items():
        for (cb,gb),vb in b.items():
            word=ga+gb
            if len(set(word))!=len(word): continue
            sign=(-1)**sum(word[j]>word[k] for j in range(len(word)) for k in range(j+1,len(word)))
            key=tuple(sorted(ca+cb)),tuple(sorted(word))
            out[key]=za(out.get(key,Z),zm(sign,zm(va,vb)))
            if not nz(out[key]): del out[key]
    return out
def times(*factors):
    out=constant(1)
    for p in factors: out=mul(out,p)
    return out


def replace(a,rule):
    """Even derivation, preserving the original ordered odd-field slots."""
    out={}
    for (comm,grass),coeff in a.items():
        factors=[(False,v) for v in comm]+[(True,v) for v in grass]
        for j,(isodd,key) in enumerate(factors):
            replaced=rule(isodd,key)
            if replaced is None: continue
            parts=[(odd(*v) if g else bos(*v)) if k!=j else replaced
                   for k,(g,v) in enumerate(factors)]
            out=add(out,scale(coeff,times(*parts)))
    return out
def substitute_odd(a,rule):
    out={}
    for (comm,grass),coeff in a.items():
        factors=[bos(*v) for v in comm]+[rule(v) for v in grass]
        out=add(out,scale(coeff,times(*factors)))
    return out


identity=((1,0),(0,1)); eps=((0,1),(-1,0))
pauli=(((0,1),(1,0)),((0,zm(-1,I)),(I,0)),((1,0),(0,-1)))
bs=(identity,)+tuple(tuple(tuple(zm(-1,e) for e in row) for row in p) for p in pauli)
T=tuple(tuple(tuple(zm(R(1,2),e) for e in row) for row in p) for p in pauli)
F=(((R(1,2),R(2)),(R(-2),R(1,3))),((R(3),R(-1)),(R(1),R(2))))
Y=(((R(2),R(-1)),(R(1,3),R(3))),((R(-1),R(2)),(R(4),R(1))))
M=((R(2),R(1)),(R(1),R(5)))
NF=2


rho=total(times(bos('pd',a),bos('p',a)) for a in range(2))
def bilinear(mu,a):
    return total(scale(zm(F[alpha][i],bs[mu][s][t]),times(odd('lb',alpha,a,s),odd('n',i,t)))
                 for alpha,i,s,t in product(range(NF),range(NF),range(2),range(2)))
b=tuple(tuple(bilinear(mu,a) for a in range(2)) for mu in range(4))
j=tuple(total(times(b[mu][a],bos('p',a)) for a in range(2)) for mu in range(4))
A=total(times(rho,b[mu][a],bos('dp',mu,a)) for mu,a in product(range(4),range(2)))
B=total(times(j[mu],bos('pd',a),bos('dp',mu,a)) for mu,a in product(range(4),range(2)))
D=total(times(j[mu],bos('dpd',mu,a),bos('p',a)) for mu,a in product(range(4),range(2)))


def partial(a,mu):
    def rule(isodd,key):
        tag=key[0]
        if not isodd and tag=='p':
            q=key[1]
            return add(bos('dp',mu,q),
                       total(times(bos('ac',mu,q,c),bos('p',c)) for c in range(2)),zm(-1,I))
        if not isodd and tag=='pd':
            q=key[1]
            return add(bos('dpd',mu,q),total(times(bos('pd',c),bos('ac',mu,c,q)) for c in range(2)),I)
        if isodd and tag=='lb':
            _,alpha,q,s=key
            return add(odd('dlb',mu,alpha,q,s),
                       total(times(odd('lb',alpha,c,s),bos('ac',mu,c,q)) for c in range(2)),I)
        if isodd and tag=='n': return odd('dn',mu,key[1],key[2])
        raise ValueError('first jet derivative outside declared inputs '+str(key))
    return replace(a,rule)


divj=total(partial(j[mu],mu) for mu in range(4))
divrhoj=total(partial(times(rho,j[mu]),mu) for mu in range(4))
def kbar(alpha,a,t):
    return total(scale(zm(-1,zm(I,bs[mu][s][t])),odd('dlb',mu,alpha,a,s))
                 for mu,s in product(range(4),range(2)))
def an(i,s):
    return total(scale(zm(I,bs[mu][s][t]),odd('dn',mu,i,t)) for mu,t in product(range(4),range(2)))
rowdiv=total(scale(F[alpha][i],times(kbar(alpha,a,t),bos('p',a),odd('n',i,t)))
             for alpha,i,a,t in product(range(NF),range(NF),range(2),range(2)))
coldiv=total(scale(F[alpha][i],times(odd('lb',alpha,a,s),bos('p',a),an(i,s)))
             for alpha,i,a,s in product(range(NF),range(NF),range(2),range(2)))
scalargrad=total(times(b[mu][a],bos('dp',mu,a)) for mu,a in product(range(4),range(2)))


def scalar_source(independent=True):
    if independent:
        sl=lambda alpha,a,t: odd('sb',alpha,a,t)
        sr=lambda i,s: odd('sr',i,s)
        nr=lambda i,s: odd('nr',i,s)
    else: sl,sr,nr=bar_sl,s_r,n_r
    first=total(scale(F[alpha][i],times(sl(alpha,a,t),bos('p',a),odd('n',i,t)))
                for alpha,i,a,t in product(range(NF),range(NF),range(2),range(2)))
    second=total(scale(F[alpha][i],times(odd('lb',alpha,a,s),bos('p',a),
                 add(total(scale(M[i][k],nr(k,s)) for k in range(NF)),sr(i,s))))
                 for alpha,i,a,s in product(range(NF),range(NF),range(2),range(2)))
    return times(rho,add(first,second,-1))


residual=times(rho,add(
    total(scale(F[alpha][i],times(odd('eb',alpha,a,t),bos('p',a),odd('n',i,t)))
          for alpha,i,a,t in product(range(NF),range(NF),range(2),range(2))),
    total(scale(F[alpha][i],times(odd('lb',alpha,a,s),bos('p',a),odd('er',i,s)))
          for alpha,i,a,s in product(range(NF),range(NF),range(2),range(2))),-1))
def eom_rule(key):
    if key[0]=='dlb' and key[1]==0:
        _,_,alpha,a,t=key
        return add(scale(I,add(odd('sb',alpha,a,t),odd('eb',alpha,a,t))),
                   total(scale(bs[mu][s][t],odd('dlb',mu,alpha,a,s))
                         for mu,s in product(range(1,4),range(2))),-1)
    if key[0]=='dn' and key[1]==0:
        _,_,i,s=key
        sources=total((total(scale(M[i][k],odd('nr',k,s)) for k in range(NF)),odd('sr',i,s),odd('er',i,s)))
        return add(scale(zm(-1,I),sources),total(scale(bs[mu][s][t],odd('dn',mu,i,t))
                   for mu,t in product(range(1,4),range(2))),-1)
    return odd(*key)


def literal_scalar(g2,gp2):
    hyper=total(scale(zm(I,R(1,4)),times(j[mu],add(times(bos('pd',a),bos('dp',mu,a)),
                                                   times(bos('dpd',mu,a),bos('p',a)),-1)))
                for mu,a in product(range(4),range(2)))
    weak={}
    for mu,q,a,c,d,e in product(range(4),range(3),range(2),range(2),range(2),range(2)):
        ji=scale(T[q][a][c],times(b[mu][a],bos('p',c)))
        scalar=scale(zm(I,T[q][d][e]),add(times(bos('pd',d),bos('dp',mu,e)),
                                                        times(bos('dpd',mu,d),bos('p',e)),-1))
        weak=add(weak,times(ji,scalar))
    return scale(-R(1,6),add(scale(gp2,hyper),scale(g2,weak)))


counts={}; mutations={}
def check(category,label,a,c):
    if a!=c: raise RuntimeError('FAILED '+category+': '+label)
    counts[category]=counts.get(category,0)+1
def reject(label,a,c):
    if a==c: raise RuntimeError('MUTATION SURVIVED '+label)
    mutations[label]=mutations.get(label,0)+1


def rank(rows):
    rows=[list(map(R,row)) for row in rows]; r=0
    for column in range(len(rows[0]) if rows else 0):
        pivot=next((k for k in range(r,len(rows)) if rows[k][column]),None)
        if pivot is None: continue
        rows[r],rows[pivot]=rows[pivot],rows[r]
        divisor=rows[r][column]; rows[r]=[v/divisor for v in rows[r]]
        for k in range(len(rows)):
            if k==r or not rows[k][column]: continue
            factor=rows[k][column]; rows[k]=[v-factor*w for v,w in zip(rows[k],rows[r])]
        r+=1
    return r


def matrix_product(a,c):
    return [[sum_complex(zm(a[i][k],c[k][j]) for k in range(len(c)))
             for j in range(len(c[0]))] for i in range(len(a))]
def sum_complex(values):
    out=Z
    for v in values: out=za(out,v)
    return out


def weak_generator(a,t):
    def rule(isodd,key):
        if isodd: raise ValueError('formal weak module uses an even covector')
        tag=key[0]; q=key[-1]
        row=tag in ('bb','pd','dpd')
        prefix=key[:-1]
        return total(scale(-t[c][q] if row else t[q][c],bos(*(prefix+(c,)))) for c in range(2))
    return replace(a,rule)


def weak_kernel(conjugate):
    if not conjugate:
        basis=[times(bos('bb',a),bos('p',c),bos('pd',d),bos('dp',0,e))
               for a,c,d,e in product(range(2),repeat=4)]
    else:
        basis=[times(bos('bb',a),bos('p',c),bos('p',d),bos('dpd',0,e))
               for a,(c,d),e in product(range(2),combinations_with_replacement(range(2),2),range(2))]
    keys=[next(iter(p)) for p in basis]; row_index={k:j for j,k in enumerate(keys)}
    constraints=[]
    generators=(((0,1),(0,0)),((0,0),(1,0)),((1,0),(0,-1)))
    for t in generators:
        matrix=[[R(0) for _ in basis] for _ in basis]
        for column,p in enumerate(basis):
            for key,v in weak_generator(p,t).items():
                if v[1]: raise RuntimeError('complex weak Lie coefficient')
                matrix[row_index[key]][column]=v[0]
        constraints.extend(matrix)
    r=rank(constraints)
    return {'derivative_on':'phi_dagger' if conjugate else 'phi','unknowns':len(basis),
            'constraint_rank':r,'invariant_dimension':len(basis)-r}


def star_variable(key,isodd):
    mapping={'p':'pd','pd':'p','dp':'dpd','dpd':'dp','l':'lb','lb':'l','n':'nb','nb':'n'}
    if key[0] not in mapping: raise ValueError('undeclared conjugate variable '+str(key))
    return (mapping[key[0]],)+key[1:]
def dagger(a,reverse=True):
    out={}
    for (comm,grass),coeff in a.items():
        word=grass[::-1] if reverse else grass
        factors=[bos(*star_variable(v,False)) for v in comm]+[odd(*star_variable(v,True)) for v in word]
        out=add(out,scale(zs(coeff),times(*factors)))
    return out


def n_r(i,s): return total(scale(-eps[s][t],odd('nb',i,t)) for t in range(2))
def bar_nr(i,t): return total(scale(eps[s][t],odd('n',i,s)) for s in range(2))
def bar_x(alpha,s): return total(times(odd('lb',alpha,a,s),bos('p',a)) for a in range(2))
def x_c(alpha,s): return total(scale(-eps[s][t],bar_x(alpha,t)) for t in range(2))
def bar_sl(alpha,a,t):
    return total(scale(zs(Y[alpha][k]),times(bar_nr(k,t),bos('pd',a))) for k in range(NF))
def s_r(i,s): return total(scale(Y[beta][i],x_c(beta,s)) for beta in range(NF))
def sn(j,i): return total(times(bar_nr(j,t),odd('n',i,t)) for t in range(2))
def sx(alpha,beta): return total(times(bar_x(alpha,s),x_c(beta,s)) for s in range(2))


def compact_sources(wrong_transpose=False):
    sterile={}; active={}; mass={}
    for j,i in product(range(NF),repeat=2):
        coeff=z(0)
        for alpha in range(NF): coeff=za(coeff,zm(zs(Y[alpha][j]),F[alpha][i]))
        sterile=add(sterile,scale(coeff,sn(j,i)))
    for alpha,beta in product(range(NF),repeat=2):
        coeff=z(0)
        for i in range(NF):
            coeff=za(coeff,zm(F[alpha][i],zs(Y[i][beta]) if wrong_transpose else Y[beta][i]))
        active=add(active,scale(coeff,sx(alpha,beta)))
    for alpha,j in product(range(NF),repeat=2):
        coeff=z(0)
        for i in range(NF): coeff=za(coeff,zm(F[alpha][i],M[i][j]))
        mass=add(mass,scale(coeff,total(times(bar_x(alpha,s),n_r(j,s)) for s in range(2))))
    return add(times(rho,rho,sterile),times(rho,add(mass,active)),-1)


def run():
    gamma0=[[z(1 if i%2==k%2 and i//2!=k//2 else 0) for k in range(4)] for i in range(4)]
    gamma2=[[z(0) for _ in range(4)] for _ in range(4)]
    for s,t in product(range(2),repeat=2):
        gamma2[s][t+2]=z(pauli[1][s][t]); gamma2[s+2][t]=zm(-1,pauli[1][s][t])
    cc=[[zm(I,v) for v in row] for row in matrix_product(gamma2,gamma0)]
    for s,t in product(range(4),repeat=2):
        expected=eps[s][t] if s<2 and t<2 else -eps[s-2][t-2] if s>=2 and t>=2 else 0
        check('charge_conjugation_matrix','C=i gamma2 gamma0',cc[s][t],z(expected))
    kernel_phi=weak_kernel(False); kernel_dagger=weak_kernel(True)
    check('weak_invariant_kernels','Dphi dimension',kernel_phi['invariant_dimension'],2)
    check('weak_invariant_kernels','Dphi dagger dimension',kernel_dagger['invariant_dimension'],1)
    # Explicit actual spinor/Grassmann component rows, not just weak tensors.
    oddword=tuple(sorted((('lb',0,0,0),('n',0,0))))
    witness_keys=[(tuple(sorted((('pd',1),('p',1),('dp',0,0)))),oddword),
                  (tuple(sorted((('pd',1),('p',0),('dp',0,1)))),oddword),
                  (tuple(sorted((('p',0),('p',1),('dpd',0,1)))),oddword)]
    component_rows=[]
    for row,key in enumerate(witness_keys):
        entries=[]
        for column,p in enumerate((A,B,D)):
            value=p.get(key,Z); check('actual_component_basis_minor',str((row,column)),value,F[0][0] if row==column else Z)
            entries.append([str(v) for v in value])
        component_rows.append(entries)
    for t in (((0,1),(0,0)),((0,0),(1,0)),((1,0),(0,-1))):
        aw=total(times(rho,bos('bb',a),bos('dp',0,a)) for a in range(2))
        bw=total(times(bos('bb',a),bos('p',a),bos('pd',c),bos('dp',0,c)) for a,c in product(range(2),repeat=2))
        dw=total(times(bos('bb',a),bos('p',a),bos('dpd',0,c),bos('p',c)) for a,c in product(range(2),repeat=2))
        for label,p in (('A',aw),('B',bw),('D',dw)): check('weak_invariant_representatives',label,weak_generator(p,t),{})
    check('complete_neutral_current_divergence','row scalar and sterile derivatives',divj,add(scalargrad,scale(I,add(rowdiv,coldiv,-1))))
    check('rho_current_product_rule','complete neutral product',divrhoj,add(times(rho,divj),add(B,D)))
    reject('discard_scalar_derivative_inside_divergence',divj,scale(I,add(rowdiv,coldiv,-1)))
    reject('replace_complete_divergence_by_sterile_derivative',divj,scale(zm(-1,I),coldiv))
    for g2,gp2 in ((R(1),R(0)),(R(0),R(1)),(R(5,3),R(7,5))):
        c=(g2+gp2)/24
        literal=literal_scalar(g2,gp2)
        frame=scale(zm(I,R(1,24)),total((scale(-2*g2,A),scale(g2-gp2,B),scale(g2+gp2,D))))
        check('literal_scalar_current_reconstruction','Pauli plus hypercharge',literal,frame)
        deriv=total((scale(zm(-1,zm(I,(3*g2+gp2)/24)),A),scale(zm(-1,zm(I,gp2/12)),B)))
        before=add(deriv,scale(c,times(rho,add(rowdiv,coldiv,-1))))
        check('off_shell_IBP_identity','retained total derivative',literal,add(before,scale(zm(I,c),divrhoj)))
        red=add(deriv,scale(c,scalar_source()))
        transported=add(red,scale(c,residual))
        check('sourced_equation_residual','independent nonzero residuals',substitute_odd(before,eom_rule),transported)
        check('full_off_shell_certificate','all jets and source residuals',substitute_odd(literal,eom_rule),
              add(transported,scale(zm(I,c),substitute_odd(divrhoj,eom_rule))))
        reject('omit_total_divergence',literal,before)
        reject('omit_lepton_and_Majorana_sources',transported,add(deriv,scale(c,residual)))
        sterile_error=times(rho,total(scale(F[alpha][i],times(odd('lb',alpha,a,s),bos('p',a),odd('er',i,s)))
                   for alpha,i,a,s in product(range(NF),range(NF),range(2),range(2))))
        wrong_residual=add(residual,scale(2,sterile_error))
        reject('reverse_Majorana_equation_sign',substitute_odd(before,eom_rule),add(red,scale(c,wrong_residual)))
        # The residual-cancelling local map has no scalar or gauge shift.
        ell=total(scale(zm(c,F[alpha][i]),times(rho,odd('eb',alpha,a,t),bos('p',a),odd('n',i,t)))
                  for alpha,i,a,t in product(range(NF),range(NF),range(2),range(2)))
        barxer=total(scale(zm(-c,F[alpha][i]),times(rho,odd('lb',alpha,a,s),bos('p',a),odd('er',i,s)))
                    for alpha,i,a,s in product(range(NF),range(NF),range(2),range(2)))
        check('additional_field_map_rows','ell and bar x',add(ell,barxer),scale(c,residual))
        reject('extra_Majorana_map_half',add(ell,scale(R(1,2),barxer)),scale(c,residual))
    for i,jj in product(range(NF),repeat=2):
        check('Grassmann_scalar_slot_symmetry','sterile',sn(jj,i),sn(i,jj))
        expected=total(times(odd('nb',i,s),n_r(jj,s)) for s in range(2))
        check('Grassmann_conjugate_order','whole sterile bilinear',dagger(sn(jj,i)),expected)
        reject('conjugate_without_reversing_odd_slots',dagger(sn(jj,i)),dagger(sn(jj,i),False))
    for alpha,beta in product(range(NF),repeat=2):
        check('Grassmann_scalar_slot_symmetry','active',sx(alpha,beta),sx(beta,alpha))
    check('explicit_Yukawa_source_orientation','Ydag F, FM, FYtranspose',scalar_source(False),compact_sources())
    reject('active_source_uses_Ydag_instead_of_Ytranspose',scalar_source(False),compact_sources(True))
    # Antisymmetric coefficient components disappear in scalar, unlike tensor, slots.
    anti=total(scale((-1 if a<c else 1) if a!=c else 0,sx(a,c)) for a,c in product(range(NF),repeat=2))
    check('scalar_coefficient_projection','active antisymmetric kernel',anti,{})
    sterile_anti=total(scale((-1 if a<c else 1) if a!=c else 0,sn(a,c)) for a,c in product(range(NF),repeat=2))
    check('scalar_coefficient_projection','sterile antisymmetric kernel',sterile_anti,{})
    reject('discard_symmetric_scalar_channel',sx(0,0),{})
    for i,t in product(range(NF),range(2)):
        x=scale(-1,total(scale(zs(F[alpha][i]),times(rho,bos('pd',a),odd('l',alpha,a,t)))
                         for alpha,a in product(range(NF),range(2))))
        barx=scale(-1,total(scale(F[alpha][i],times(rho,odd('lb',alpha,a,t),bos('p',a)))
                            for alpha,a in product(range(NF),range(2))))
        check('additional_field_map_adjoint','x versus bar x',dagger(x),barx)
        x_components=[scale(-1,total(scale(zs(F[alpha][i]),times(rho,bos('pd',a),odd('l',alpha,a,u)))
                          for alpha,a in product(range(NF),range(2)))) for u in range(2)]
        right=[total(scale(-eps[s][u],dagger(x_components[u])) for u in range(2)) for s in range(2)]
        cc_left=total(scale(eps[t][s],dagger(right[s])) for s in range(2))
        check('additional_field_map_Majorana_completion','charge conjugate has original left slot',cc_left,x)
    result={'status':'PASS','method':'exact rational complex weak-invariant kernels, commuting scalar jets and exterior Weyl components',
            'tests_by_category':counts,'positive_checks':sum(counts.values()),
            'mutation_rejections':mutations,'negative_checks':sum(mutations.values()),
            'total_checks':sum(counts.values())+sum(mutations.values()),
            'weak_kernels':[kernel_phi,kernel_dagger],
            'component_minor':{'basis':['A=rho b Dphi','B=j phi_dagger Dphi','D=j (Dphi)_dagger phi'],
                               'matrix_real_imaginary_fraction_strings':component_rows,
                               'determinant':'F_(00)^3; F_(00)=1/2+2i, nonzero'},
            'flavor_components':{'active':NF,'sterile':NF},'gauge_coefficient_trials':3,
            'scope':['only scalar-current contribution of selected sourced packet',
                     'full off-shell total derivative plus sourced lepton/Majorana residual certificate',
                     'four-dimensional fixed chirality; no physical total coefficient or beta function',
                     'basis independence is in this scalar-jet sector, not full EFT EOM closure']}
    target=Path(__file__).with_name('INDEPENDENT_SCALAR_AUDIT_RESULTS.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': run()
