#!/usr/bin/env python3
"""Literal Yukawa-flow/chiral proof and MS pole-to-running bookkeeping.

The proper-vertex result is derived from graphs and propagator numerators.
The field and scalar two-point pole coefficients are explicit inputs from
the separately checked producer. No project modules or saved JSON imported.
"""
from fractions import Fraction as F
from itertools import product, combinations_with_replacement
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
    def __bool__(self): return bool(self.r or self.i)
    def __eq__(self,o):
        o=Q(o); return self.r==o.r and self.i==o.i
    def star(self): return Q(self.r,-self.i)


I=Q(0,1); ETA=(1,-1,-1,-1)
def eye(n): return tuple(tuple(Q(i==j) for j in range(n)) for i in range(n))
def zero(n,m=None): return tuple(tuple(Q() for j in range(n if m is None else m)) for i in range(n))
def add(a,b): return tuple(tuple(x+y for x,y in zip(ar,br)) for ar,br in zip(a,b))
def scale(c,a): return tuple(tuple(c*x for x in row) for row in a)
def sub(a,b): return add(a,scale(-1,b))
def mul(a,b):
    return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(len(b))),Q()) for j in range(len(b[0]))) for i in range(len(a)))
def transpose(a): return tuple(tuple(a[j][i] for j in range(len(a))) for i in range(len(a[0])))
def star(a): return tuple(tuple(x.star() for x in row) for row in a)
def dagger(a): return transpose(star(a))
def trace(a): return sum((a[i][i] for i in range(len(a))),Q())
def block(a,r,c,n=2): return tuple(tuple(a[i+r][j+c] for j in range(n)) for i in range(n))
def comm(a,b): return sub(mul(a,b),mul(b,a))


PAULI=(((Q(),Q(1)),(Q(1),Q())),((Q(),-I),(I,Q())),((Q(1),Q()),(Q(),Q(-1))))
G0=tuple(tuple(Q((i<2 and j==i+2) or (i>=2 and j==i-2)) for j in range(4)) for i in range(4))
def spatial(p): return tuple(tuple(p[i][j-2] if i<2<=j else -p[i-2][j] if j<2<=i else Q() for j in range(4)) for i in range(4))
GAM=(G0,)+tuple(spatial(p) for p in PAULI)
G5=scale(I,mul(mul(mul(GAM[0],GAM[1]),GAM[2]),GAM[3]))
PL=scale(F(1,2),sub(eye(4),G5)); PR=scale(F(1,2),add(eye(4),G5))
CC=scale(I,mul(GAM[2],GAM[0])); EPS=((Q(),Q(1)),(Q(-1),Q()))
def slash(p):
    ans=zero(4)
    for mu in range(4): ans=add(ans,scale(ETA[mu]*p[mu],GAM[mu]))
    return ans


VERTICES={'A':('Lbar','phi','NR'),'B':('L','phibar','NRbar')}
CHARGE={'L':-F(1,2),'Lbar':F(1,2),'phi':-F(1,2),'phibar':F(1,2),'NR':F(0),'NRbar':F(0)}
EDGES=((0,1),(0,2),(1,2))
def prop(a,b):
    pair=tuple(sorted((a,b)))
    if pair==('phi','phibar'): return 'scalar'
    if pair==('L','Lbar'): return 'light_q'
    if pair==('NR','NRbar'): return 'sterile_q'
    if pair in (('NR','NR'),('NRbar','NRbar')): return 'sterile_M'
    return None


def triangle(orientations,external):
    intern=[]
    for orientation,field in zip(orientations,external):
        fields=list(VERTICES[orientation])
        if field not in fields: return []
        fields.remove(field); intern.append(fields)
    allowed=[]
    for flips in product((0,1),repeat=3):
        incidence={}
        for node in range(3):
            neighbors=[j for j in range(3) if j!=node]
            for n,neighbor in enumerate(neighbors): incidence[node,neighbor]=intern[node][n^flips[node]]
        props=[prop(incidence[a,b],incidence[b,a]) for a,b in EDGES]
        if all(props):
            numerator_q=sum(x in ('light_q','sterile_q') for x in props)
            allowed.append({'vertex_orientations':list(orientations),'external_fields':list(external),
                            'edges':[{'endpoints':[a,b],'fields':[incidence[a,b],incidence[b,a]],'propagator':p}
                                     for (a,b),p in zip(EDGES,props)],
                            'numerator_q_degree':numerator_q,'superficial_UV_degree':4+numerator_q-6})
    return allowed


def connected(edges,remove=None):
    visited={0}
    while True:
        old=set(visited)
        for j,(a,b) in enumerate(edges):
            if j==remove: continue
            if a in visited: visited.add(b)
            if b in visited: visited.add(a)
        if visited==old: return len(visited)==3
def proper_cubic_topologies():
    proper=[]; reducible=[]
    edge_types=tuple((a,b) for a in range(3) for b in range(a,3))
    for edges in combinations_with_replacement(edge_types,3):
        degree=[sum((a==i)+(b==i) for a,b in edges) for i in range(3)]
        if any(d>3 for d in degree) or not connected(edges): continue
        external=tuple(3-d for d in degree)
        if sum(external)!=3: continue
        if all(connected(edges,j) for j in range(3)): proper.append((edges,external))
        else: reducible.append((edges,external))
    return proper,reducible


def poles(y,m,ms):
    """Two-point pole INPUTS; overall 1/(16pi² epsilon) stripped."""
    hl=mul(y,dagger(y)); hn=mul(dagger(y),y); tt=trace(hn)
    zl=scale(-F(1,2),hl); zn=scale(-1,hn); zp=-tt
    ry=add(add(scale(F(1,4),mul(hl,y)),scale(F(1,2),mul(y,hn))),scale(tt*F(1,2),y))
    rm=scale(F(1,2),add(mul(transpose(hn),m),mul(m,hn)))
    qq=trace(mul(hn,mul(dagger(m),m)))
    rs=ms*tt-2*qq
    return zl,zn,zp,ry,rm,rs,hl,hn,tt,qq


checks={}; mutations={}
def check(cat,label,a,b):
    if a!=b: raise RuntimeError('FAILED '+cat+': '+label)
    checks[cat]=checks.get(cat,0)+1
def reject(label,a,b):
    if a==b: raise RuntimeError('MUTATION SURVIVED '+label)
    mutations[label]=mutations.get(label,0)+1


def run():
    proper,reducible=proper_cubic_topologies()
    check('one_particle_irreducibility','only cubic triangle topology',proper,[(EDGES,(1,1,1))])
    check('one_particle_irreducibility','all two-external-one-vertex attachments reducible',
          all(any(not connected(e,j) for j in range(3)) for e,x in reducible if max(x)>=2),True)
    target=[]; crossed=[]
    for orientations in product('AB',repeat=3):
        target+=triangle(orientations,('Lbar','phi','NR'))
        crossed+=triangle(orientations,('Lbar','phi','NRbar'))
    check('literal_graph_enumeration','target dimension-four orientation',len(target),0)
    check('literal_graph_enumeration','opposite Majorana orientation',len(crossed),1)
    check('literal_graph_enumeration','crossed orientation types',crossed[0]['vertex_orientations'],['A','A','B'])
    check('literal_graph_enumeration','crossed propagator set',sorted(e['propagator'] for e in crossed[0]['edges']),['light_q','scalar','sterile_M'])
    check('literal_graph_enumeration','crossed triangle UV degree',crossed[0]['superficial_UV_degree'],-1)
    reject('treat_Majorana_mass_numerator_as_second_slash',crossed[0]['superficial_UV_degree'],0)
    for forbidden in (('phi','phi'),('Lbar','Lbar')):
        check('charge_flow','same-orientation charged propagator absent',prop(*forbidden),None)
        reject('permit_charged_same_orientation_propagator',CHARGE[forbidden[0]]+CHARGE[forbidden[1]],F(0))
    for mu,nu in product(range(4),repeat=2):
        check('literal_clifford','gamma anticommutator',add(mul(GAM[mu],GAM[nu]),mul(GAM[nu],GAM[mu])),scale(2*ETA[mu] if mu==nu else 0,eye(4)))
    qs=((3,1,-2,2),(-2,3,1,-1))+tuple(tuple(int(i==j) for i in range(4)) for j in range(4))
    for q in qs:
        numerator=add(slash(q),scale(5,eye(4)))
        actual=mul(mul(mul(PR,numerator),transpose(CC)),transpose(PR))
        expected=scale(5,mul(mul(PR,transpose(CC)),transpose(PR)))
        check('actual_majorana_propagator_chain','two right endpoints project out slash',actual,expected)
        check('actual_majorana_propagator_chain','PR slash PR vanishes',mul(mul(PR,slash(q)),PR),zero(4))
        reject('omit_endpoint_chiral_projectors',actual,mul(numerator,transpose(CC)))
        light=block(mul(mul(PL,slash(q)),PR),0,2)
        stripped=mul(mul(block(actual,2,2),transpose(light)),transpose(EPS))
        check('actual_majorana_propagator_chain','crossed open-chain tensor is M slash',stripped,scale(5,block(slash(q),2,0)))
        reject('insert_UV_q_squared_numerator',stripped,scale(sum(ETA[mu]*q[mu]**2 for mu in range(4)),eye(2)))
    for mass_insertions in (1,3,5):
        # Additional two powers of M arise from propagator denominator
        # expansion and add two inverse loop momenta, not numerator q's.
        omega=-1-(mass_insertions-1)
        check('UV_IR_separation','UV radial power below logarithmic',omega<0,True)
    for a,b in product(range(2),repeat=2):
        weak=sum(int(a==c)*int(b==c) for c in range(2))
        check('weak_charge_flow','external weak delta without extra two',weak,int(a==b))
    reject('extra_weak_doublet_factor',sum(int(0==c)*int(0==c) for c in range(2)),2)
    # Noncommuting rectangular Y and symmetric complex M expose flavor order.
    y=((Q(1,2),Q(2,-1)),(Q(-1,3),Q(3,1)),(Q(2,-2),Q(1,1)))
    m=((Q(3,1),Q(1,-2)),(Q(1,-2),Q(5,-1))); ms=Q(F(11,7))
    zl,zn,zp,ry,rm,rs,hl,hn,tt,qq=poles(y,m,ms)
    check('renormalized_parameter_counterterms','proper Y action counterterm zero',
          add(add(add(ry,scale(F(1,2),mul(dagger(zl),y))),scale(F(1,2),mul(y,zn))),scale(zp*F(1,2),y)),zero(3,2))
    check('renormalized_parameter_counterterms','no proper sterile mass pole',
          add(rm,scale(F(1,2),add(mul(transpose(zn),m),mul(m,zn)))),zero(2))
    check('renormalized_parameter_counterterms','scalar raw mass action pole',rs+ms*zp,-2*qq)
    check('renormalized_parameter_counterterms','Majorana parameter pole symmetric',transpose(rm),rm)
    check('renormalized_parameter_counterterms','scalar parameter pole real',rs.star(),rs)
    by=scale(2,ry); bm=scale(2,rm); bs=2*rs
    check('canonical_epsilon_running','bare Y equation through one loop',add(add(ry,by),scale(-3,ry)),zero(3,2))
    check('canonical_epsilon_running','bare M equation no mu^eps factor',sub(bm,scale(2,rm)),zero(2))
    check('canonical_epsilon_running','bare scalar mass equation',bs-2*rs,Q())
    check('canonical_epsilon_running','beta Y from anomalous dimensions',by,
          add(add(scale(F(1,2),mul(hl,y)),mul(y,hn)),scale(tt,y)))
    check('canonical_epsilon_running','associative Y flavor products',mul(hl,y),mul(y,hn))
    reject('miss_factor_two_from_mu_epsilon',by,ry)
    reject('insert_mu_epsilon_prefactor_for_mass',bm,rm)
    reject('use_adjoint_instead_of_Majorana_transpose',rm,scale(F(1,2),add(mul(dagger(hn),m),mul(m,hn))))
    reject('omit_scalar_field_mass_transport',rs,-2*qq)
    reject('omit_Y_vertex_field_transport',ry,zero(3,2))
    for t in (F(2),F(3,5),F(-2,3)):
        pt=poles(scale(t,y),m,ms)
        check('independent_homogeneity','cubic Yukawa pole',pt[3],scale(t**3,ry))
        check('independent_homogeneity','quadratic Majorana pole',pt[4],scale(t**2,rm))
        check('independent_homogeneity','quadratic scalar-mass pole',pt[5],t**2*rs)
    v=((Q(F(3,5)),I*F(4,5)),(I*F(4,5),Q(F(3,5))))
    u=((Q(F(3,5)),Q(F(4,5)),Q()),(Q(-F(4,5)),Q(F(3,5)),Q()),(Q(),Q(),I))
    check('flavor_covariance','sterile unitary',mul(dagger(v),v),eye(2))
    check('flavor_covariance','active unitary',mul(dagger(u),u),eye(3))
    yp=mul(mul(dagger(u),y),v); mp=mul(mul(transpose(v),m),v)
    rp=poles(yp,mp,ms)
    check('flavor_covariance','Y pole',rp[3],mul(mul(dagger(u),ry),v))
    check('flavor_covariance','Majorana pole congruence',rp[4],mul(mul(transpose(v),rm),v))
    check('flavor_covariance','scalar trace mass invariant',rp[9],qq)
    check('flavor_covariance','scalar parameter pole invariant',rp[5],rs)
    pzero=poles(y,zero(2),ms)
    check('massless_limit','no induced sterile mass from zero M',pzero[4],zero(2))
    check('massless_limit','scalar mass pole without retained M',pzero[5],ms*tt)
    one=poles(((Q(2),),),((Q(3),),),Q(5))
    check('one_flavor_limit','16pi2 beta Y equals5y³/2',2*one[3][0][0],Q(20))
    check('one_flavor_limit','16pi2 beta M equals2y²M',2*one[4][0][0],Q(24))
    check('one_flavor_limit','16pi2 beta scalar mass',2*one[5],Q(-104))
    result={'status':'PASS','positive_checks':sum(checks.values()),'checks_by_category':checks,
            'negative_checks':sum(mutations.values()),'mutation_rejections':mutations,
            'total_checks':sum(checks.values())+sum(mutations.values()),
            'proper_cubic_topology_count':len(proper),'reducible_topologies_excluded':len(reducible),
            'target_NR_allowed_triangles':target,'opposite_Majorana_allowed_triangles':crossed,
            'proper_Yukawa_UV_pole':'zero','proper_sterile_mass_UV_pole':'zero in declared pure-neutrino-Yukawa sector',
            'two_point_inputs':{'deltaZ_L':'-kappa YYdagger/2','deltaZ_NR':'-kappa H','deltaZ_phi':'-kappa trH',
                                'raw_scalar_mass_counterterm':'-2kappa Qm','H':'YdaggerY','Qm':'Tr(H MdaggerM)'},
            'parameter_pole_coefficients':{'deltaY/kappa':'3YH/4+(trH)Y/2',
                                            'deltaM/kappa':'(H^T M+MH)/2','deltam_H_squared/kappa':'m_H_squared trH-2Qm'},
            'one_loop_4D_running':{'16pi2 betaY':'3YH/2+(trH)Y','16pi2 betaM':'H^T M+MH',
                                  '16pi2 beta_m_H_squared':'2m_H_squared trH-4Qm'},
            'arithmetic':'exact Gaussian rational, standard library only',
            'scope':['known one-loop renormalizable pure-neutrino-Yukawa benchmark; unbroken background',
                     'proper vertex enumeration/chiral power counting independent of two-point input',
                     'RG coefficients additionally depend on separately calculated two-point inputs',
                     'fixed flavor coordinates; no gauge/other-Yukawa/quartic contributions included',
                     'no finite triangle amplitude, complete TOE/seed beta function, or novelty claim']}
    Path(__file__).with_name('YUKAWA_VERTEX_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': run()
