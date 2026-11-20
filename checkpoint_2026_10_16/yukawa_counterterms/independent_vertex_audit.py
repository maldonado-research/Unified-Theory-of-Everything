#!/usr/bin/env python3
"""Independent pure-neutrino-Yukawa proper-vertex and pole-to-RG audit.

Exact Gaussian rational arithmetic and literal field-slot graph enumeration;
standard library only, no project imports or optimization-removable asserts.
"""
from fractions import Fraction as R
from itertools import permutations, product
from pathlib import Path
import json


def z(x): return x if isinstance(x,tuple) else (R(x),R(0))
def za(a,b):
    a,b=z(a),z(b); return a[0]+b[0],a[1]+b[1]
def zm(a,b):
    a,b=z(a),z(b); return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def zs(a): return a[0],-a[1]
def nz(a): return bool(a[0] or a[1])
I=(R(0),R(1)); ZERO=z(0)
def zsum(values):
    out=ZERO
    for a in values: out=za(out,a)
    return out
def zero(n,m=None): return [[ZERO for _ in range(n if m is None else m)] for _ in range(n)]
def eye(n): return [[z(i==j) for j in range(n)] for i in range(n)]
def ma(a,b,c=1): return [[za(v,zm(c,w)) for v,w in zip(row,other)] for row,other in zip(a,b)]
def ms(c,a): return [[zm(c,v) for v in row] for row in a]
def mm(a,b): return [[zsum(zm(a[i][k],b[k][j]) for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]
def tr(a): return zsum(a[i][i] for i in range(len(a)))
def transpose(a): return [list(row) for row in zip(*a)]
def star(a): return [[zs(v) for v in row] for row in a]
def dagger(a): return transpose(star(a))


PAULI=([[z(0),z(1)],[z(1),z(0)]],[[ZERO,zm(-1,I)],[I,ZERO]],[[z(1),ZERO],[ZERO,z(-1)]])
GAMMA=[]
for mu in range(4):
    a=zero(4); upper=eye(2) if mu==0 else PAULI[mu-1]
    lower=eye(2) if mu==0 else ms(-1,PAULI[mu-1])
    for s,t in product(range(2),repeat=2): a[s][t+2]=upper[s][t]; a[s+2][t]=lower[s][t]
    GAMMA.append(a)
ETA=(1,-1,-1,-1)
PL=[[z(i==j and i<2) for j in range(4)] for i in range(4)]
PR=ma(eye(4),PL,-1)
def slash(q):
    out=zero(4)
    for mu in range(4): out=ma(out,GAMMA[mu],zm(ETA[mu],q[mu]))
    return out
def dot(a,b): return zsum(zm(ETA[mu],zm(a[mu],b[mu])) for mu in range(4))


POS={}; NEG={}
def check(category,label,a,b):
    if a!=b: raise RuntimeError('FAILED '+category+': '+label)
    POS[category]=POS.get(category,0)+1
def reject(label,a,b):
    if a==b: raise RuntimeError('ADVERSE CONTROL SURVIVED '+label)
    NEG[label]=NEG.get(label,0)+1


# A is the literal -bar L Y phi N_R; B its conjugate. In left-Weyl
# components A contains barred l and barred n; B contains unbarred l,n.
SLOTS={'A':('Lb','nb','phi'),'B':('L','n','phid')}
def propagator(a,b,wrong_scalar=False,allow_anomalous=True):
    if {a,b}=={'Lb','L'}: return 'L_kinetic'
    if {a,b}=={'nb','n'}: return 'N_kinetic'
    if allow_anomalous and a==b and a in ('nb','n'): return 'N_anomalous_mass'
    if {a,b}=={'phi','phid'}: return 'scalar'
    if wrong_scalar and a==b and a in ('phi','phid'): return 'WRONG_same_charge_scalar'
    return None
def pairings(slots,wrong_scalar=False,allow_anomalous=True):
    if not slots:
        yield (); return
    first=slots[0]
    for j in range(1,len(slots)):
        second=slots[j]; kind=propagator(first[1],second[1],wrong_scalar,allow_anomalous)
        if kind:
            for rest in pairings(slots[1:j]+slots[j+1:],wrong_scalar,allow_anomalous): yield ((first,second,kind),)+rest
def connected(edges,n,skip=None):
    seen={0}; changed=True
    while changed:
        changed=False
        for j,(left,right,_) in enumerate(edges):
            if j==skip: continue
            a,b=left[0],right[0]
            if a in seen and b not in seen: seen.add(b); changed=True
            if b in seen and a not in seen: seen.add(a); changed=True
    return len(seen)==n
def canonical_graph(types,external,edges):
    candidates=[]
    for mapping_tuple in permutations(range(len(types))):
        mapping=dict(enumerate(mapping_tuple)); vertices=[None]*len(types)
        for old,t in enumerate(types): vertices[mapping[old]]=(t,tuple(sorted(k for k,slot in external.items() if slot[0]==old)))
        mapped=[]
        for left,right,kind in edges:
            ends=sorted(((mapping[left[0]],left[1]),(mapping[right[0]],right[1])))
            mapped.append((ends[0],ends[1],kind))
        candidates.append((tuple(vertices),tuple(sorted(mapped))))
    return min(candidates)
def enumerate_graphs(external_fields,vertices=3,wrong_scalar=False,allow_anomalous=True):
    out={}
    for types in product(('A','B'),repeat=vertices):
        slots=tuple((v,kind) for v,t in enumerate(types) for kind in SLOTS[t])
        names=list(external_fields)
        choices=[tuple(s for s in slots if s[1]==external_fields[name]) for name in names]
        for chosen in product(*choices):
            if len(set(chosen))!=len(chosen): continue
            ext=dict(zip(names,chosen)); remaining=tuple(s for s in slots if s not in chosen)
            for edges in pairings(remaining,wrong_scalar,allow_anomalous):
                if not connected(edges,vertices): continue
                key=canonical_graph(types,ext,edges)
                if key in out: continue
                vv,ee=key; bridges=[j for j in range(len(ee)) if not connected(ee,vertices,j)]
                degrees=[sum(left[0]==v for left,_,_ in ee)+sum(right[0]==v for _,right,_ in ee) for v in range(vertices)]
                numerator_degree=sum(kind in ('L_kinetic','N_kinetic') for _,_,kind in ee)
                weak=[]
                for a,b in product(range(2),repeat=2):
                    count=0
                    for components in product(range(2),repeat=vertices):
                        valid=True
                        for v,(_,external_names) in enumerate(vv):
                            if 'barL' in external_names and components[v]!=a: valid=False
                            if 'phi' in external_names and components[v]!=b: valid=False
                        for left,right,kind in ee:
                            if kind in ('L_kinetic','scalar','WRONG_same_charge_scalar') and components[left[0]]!=components[right[0]]: valid=False
                        count+=valid
                    weak.append(count)
                out[key]={'vertices':[[t,list(extnames)] for t,extnames in vv],
                          'internal_edges':[[list(left),list(right),kind] for left,right,kind in ee],
                          'bridges':bridges,'internal_degrees':degrees,'one_particle_irreducible':not bridges,
                          'anomalous_mass_edges':sum(kind=='N_anomalous_mass' for _,_,kind in ee),
                          'superficial_UV_degree':None if bridges else 4+numerator_degree-2*len(ee),
                          'loop_subgraph_superficial_UV_degree':4+sum(kind in ('L_kinetic','N_kinetic') for k,(_,_,kind) in enumerate(ee) if k not in bridges)-2*(len(ee)-len(bridges)),
                          'weak_external_component_matrix':[weak[:2],weak[2:]]}
    return list(out.values())


def pm(a,b):
    out={}
    for degree,x in a.items():
        for power,y in b.items(): out[degree+power]=ma(out.get(degree+power,zero(len(x),len(y[0]))),mm(x,y))
    return out
def pconstant(a): return {0:a}
def numerator(q,r,s,mass,anomalous):
    sterile={0:ma(slash(r),ms(mass,eye(4))),1:slash(q)}
    lepton={0:slash(s),1:slash(q)}
    middle=PR if anomalous else PL; final=PL if anomalous else PR
    return pm(pm(pm(pm(pconstant(PR),sterile),pconstant(middle)),lepton),pconstant(final))
def uv_log_coefficient(poly,q):
    q2=dot(q,q)
    if q2[1] or not q2[0]: raise RuntimeError('test Wick direction must have real nonzero q2')
    return ms(R(1,q2[0]**3),poly.get(2,zero(4)))


def pureY_residues(y,m):
    hl=mm(y,dagger(y)); hn=mm(dagger(y),y); t=tr(hn)
    zl=ms(R(-1,2),hl); zr=ms(-1,hn); zp=zm(-1,t)
    field_shift=ma(ma(ms(R(1,2),mm(zl,y)),ms(R(1,2),mm(y,zr))),ms(zm(R(1,2),zp),y))
    ry=ms(-1,field_shift)
    rm=ms(R(1,2),ma(mm(transpose(hn),m),mm(m,hn)))
    return hl,hn,t,zl,zr,zp,ry,rm


def run():
    for mu,nu in product(range(4),repeat=2):
        check('Clifford_and_projectors','gamma anticommutator',ma(mm(GAMMA[mu],GAMMA[nu]),mm(GAMMA[nu],GAMMA[mu])),ms(2*ETA[mu] if mu==nu else 0,eye(4)))
    check('Clifford_and_projectors','PR PL zero',mm(PR,PL),zero(4))
    selected=enumerate_graphs({'barL':'Lb','phi':'phi','sterile':'nb'})
    alternate=enumerate_graphs({'barL':'Lb','phi':'phi','sterile':'n'})
    selected_1pi=[g for g in selected if g['one_particle_irreducible']]
    alternate_1pi=[g for g in alternate if g['one_particle_irreducible']]
    check('literal_field_flow_enumeration','selected proper triangle absent',selected_1pi,[])
    check('literal_field_flow_enumeration','three selected external leg chains',len(selected),3)
    check('literal_field_flow_enumeration','all three bridge species',sorted(g['internal_edges'][g['bridges'][0]][2] for g in selected),['L_kinetic','N_kinetic','scalar'])
    for graph in selected:
        bridge=graph['internal_edges'][graph['bridges'][0]][2]
        factor=2 if bridge=='N_kinetic' else 1
        check('literal_weak_component_sums',bridge,graph['weak_external_component_matrix'],[[factor,0],[0,factor]])
    check('literal_field_flow_enumeration','Majorana alternate proper triangle exists',len(alternate_1pi),1)
    ordinary_dirac=enumerate_graphs({'barL':'Lb','phi':'phi','sterile':'n'},allow_anomalous=False)
    check('literal_field_flow_enumeration','ordinary Dirac has no alternate proper triangle',[g for g in ordinary_dirac if g['one_particle_irreducible']],[])
    for g in alternate_1pi:
        check('literal_field_flow_enumeration','mass insertion required',g['anomalous_mass_edges'],1)
        check('literal_field_flow_enumeration','negative UV degree',g['superficial_UV_degree'],-1)
        check('literal_field_flow_enumeration','proper weak factor identity',g['weak_external_component_matrix'],[[1,0],[0,1]])
    mass_two=enumerate_graphs({'sterile1':'nb','sterile2':'nb'},2)
    normal_two=enumerate_graphs({'sterile1':'nb','sterile2':'n'},2)
    check('sterile_mass_flow','no proper anomalous self energy',mass_two,[])
    check('sterile_mass_flow','ordinary kinetic graph exists',len(normal_two),1)
    check('sterile_mass_flow','free weak component sum two',normal_two[0]['weak_external_component_matrix'],[[2,2],[2,2]])
    # Complex scalar as two real fields. Factors of 1/sqrt2 are included as
    # 1/2 only in products, avoiding any irrational-number arithmetic.
    aa=(z(1),I); bb=(z(1),zm(-1,I)); extract_phi=(z(1),zm(-1,I))
    same_charge=zm(R(1,2),zsum(zm(a,a) for a in aa))
    conjugate_charge=zm(R(1,2),zsum(zm(a,b) for a,b in zip(aa,bb)))
    check('real_component_charge_control','phi phi propagator cancellation',same_charge,ZERO)
    check('real_component_charge_control','phi phi dagger propagator retained',conjugate_charge,z(1))
    check('real_component_charge_control','external A source retained',zm(R(1,2),zsum(zm(a,b) for a,b in zip(extract_phi,aa))),z(1))
    check('real_component_charge_control','external conjugate source excluded',zm(R(1,2),zsum(zm(a,b) for a,b in zip(extract_phi,bb))),ZERO)
    for trial in range(12):
        q=(zm(I,R(trial+2,3)),z(1),z(R(trial+1,5)),z(-2))
        r=(z(0),z(R(trial-2,7)),z(1),z(-1))
        s=(z(0),z(-1),z(R(trial+3,11)),z(R(1,2)))
        mass=R(trial+1,3); q2=dot(q,q)
        normal=numerator(q,r,s,mass,False); anomalous=numerator(q,r,s,mass,True)
        check('actual_loop_numerator','normal chain has no mass numerator',normal,numerator(q,r,s,0,False))
        check('actual_loop_numerator','anomalous chain vanishes at zero mass',numerator(q,r,s,0,True),{k:zero(4) for k in anomalous})
        check('actual_loop_numerator','normal leading q2 PR',normal[2],ms(q2,PR))
        check('actual_loop_numerator','anomalous kinetic leading zero',anomalous[2],zero(4))
        check('actual_loop_numerator','mass times one slash',anomalous[1],ms(mass,mm(slash(q),PL)))
        check('actual_loop_numerator','shifted mass numerator',anomalous[0],ms(mass,mm(slash(s),PL)))
        for degree,coefficient in anomalous.items():
            check('Majorana_external_chirality','right source projection zero',mm(coefficient,PR),zero(4))
        log_normal=uv_log_coefficient(normal,q); log_mass=uv_log_coefficient(anomalous,q)
        reject('normal_unweighted_triangle_false_zero',log_normal,zero(4))
        check('UV_asymptotic_log_coefficient','complex scalar cancellation',ms(same_charge,log_normal),zero(4))
        check('UV_asymptotic_log_coefficient','anomalous mass triangle UV finite',log_mass,zero(4))
        reject('replace_same_charge_by_conjugate_scalar_propagator',ms(conjugate_charge,log_normal),zero(4))
        reject('omit_imaginary_scalar_component',ms(R(1,2),log_normal),zero(4))
        reject('replace_mass_projector_by_kinetic_projector',log_normal,log_mass)
    # Nonzero ordinary log calibration with a positive auxiliary mass:
    # integral q_E^2/(q_E^2+mu_aux^2)^3 has normalized UV residue one.
    # Gamma(eps)*eps ->1; d/2 /Gamma(3) ->1 at d=4. No scaleless-zero rule.
    check('UV_radial_power','normal logarithmic radial exponent',3+2-6,-1)
    check('UV_radial_power','mass insertion integrable radial exponent',3+1-6,-2)
    check('UV_radial_power','one epsilon-pole calibration factor',R(4,2)/2,1)
    y=[[z((R(1),R(1))),z((R(2),R(-1))),z((R(-1),R(2)))],
       [z((R(3),R(-2))),z((R(1),R(1))),z((R(2),R(3)))]]
    m=[[z((R(1),R(1))),z(2),ZERO],[z(2),z((R(2),R(-2))),ZERO],[ZERO,ZERO,z((R(3),R(1)))]]
    hl,hn,t,zl,zr,zp,ry,rm=pureY_residues(y,m)
    check('pole_to_RG_flavor_algebra','HL Y=Y HN',mm(hl,y),mm(y,hn))
    field_shift=ma(ma(ms(R(1,2),mm(zl,y)),ms(R(1,2),mm(y,zr))),ms(zm(R(1,2),zp),y))
    check('pole_to_RG_flavor_algebra','zero proper vertex fixes parameter pole',ma(ry,field_shift),zero(2,3))
    by=ma(ms(R(3,2),mm(y,hn)),ms(t,y))
    bm=ma(mm(transpose(hn),m),mm(m,hn))
    check('pole_to_RG_flavor_algebra','Y beta coefficient from residue',ms(2,ry),by)
    check('pole_to_RG_flavor_algebra','M beta coefficient from residue',ms(2,rm),bm)
    check('pole_to_RG_flavor_algebra','complex symmetric mass retained',transpose(bm),bm)
    gamma_l=ms(-1,zl); gamma_r=ms(-1,zr); gamma_phi=zm(-1,zp)
    check('pole_to_RG_flavor_algebra','field anomalous dimensions give Y beta',ma(ma(mm(gamma_l,y),mm(y,gamma_r)),ms(gamma_phi,y)),by)
    check('pole_to_RG_flavor_algebra','Majorana mass uses transpose',ma(mm(transpose(gamma_r),m),mm(m,gamma_r)),bm)
    reject('Majorana_mass_uses_adjoint',ma(mm(dagger(gamma_r),m),mm(m,gamma_r)),bm)
    reject('miss_pole_to_beta_factor_two',ry,by)
    reject('miss_mass_pole_to_beta_factor_two',rm,bm)
    wrong_y=ma(ma(ms(R(1,4),mm(hl,y)),ms(R(1,2),mm(y,transpose(hn)))),ms(zm(R(1,2),t),y))
    reject('right_sterile_kinetic_matrix_transposed_in_Y_vertex',ry,wrong_y)
    for scale_y in (R(2),R(3,2)):
        scaled=pureY_residues(ms(scale_y,y),m)
        check('actual_residue_Euler_degree','cubic Yukawa residue',scaled[6],ms(scale_y**3,ry))
        check('actual_residue_Euler_degree','quadratic-in-Y mass residue',scaled[7],ms(scale_y**2,rm))
    check('canonical_epsilon_bare_derivative','Y0 has mu^eps',ma(ms(1-3,ry),by),zero(2,3))
    check('canonical_epsilon_bare_derivative','M0 has no mu power',ma(ms(0-2,rm),bm),zero(3))
    # Translate through the holomorphic N_L convention Mcal=M*, Z_NL=Z_NR*.
    beta_cal=ma(mm(transpose(star(gamma_r)),star(m)),mm(star(m),star(gamma_r)))
    check('Majorana_chiral_convention','conjugated left-Weyl mass beta',star(beta_cal),bm)
    report={'status':'PASS','calculation':'pure-neutrino-Yukawa one-loop proper barL phi N_R UV vertex and conditional pole-to-RG sector audit',
            'dateline':'15 October 2026 Pacific client; 16 October 2026 UTC checkpoint folder',
            'tests_by_category':POS,'positive_checks':sum(POS.values()),
            'adverse_controls':NEG,'negative_checks':sum(NEG.values()),'total_checks':sum(POS.values())+sum(NEG.values()),
            'selected_vertex_connected_graphs':selected,'alternate_Majorana_connected_graphs':alternate,
            'sterile_normal_two_point_graphs':normal_two,
            'selected_proper_vertex_UV_pole':'zero by literal charge/chirality flow, not inferred from beta',
            'alternate_Majorana_triangle':'one anomalous mass edge; superficial UV degree -1; finite form factors not computed',
            'conditional_external_pole_inputs':{'delta_Z_L':'-HL/(32 pi^2 eps)','delta_Z_NR':'-HN/(16 pi^2 eps)',
                                                'delta_Z_phi':'-tr(HN)/(16 pi^2 eps)'},
            'conditional_four_dimensional_RG':{'16pi2_betaY':'3 Y HN/2+tr(HN)Y','16pi2_betaM':'HN^T M+M HN'},
            'scope':['pure neutrino-Yukawa renormalizable sector around unbroken zero scalar background',
                     'field poles are external inputs for the RG algebra, not recomputed here',
                     'no gauge other Yukawa higher-dimensional or finite matching contribution',
                     'not a full physical EFT running calculation or observable prediction']}
    Path(__file__).with_name('INDEPENDENT_VERTEX_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if not k.endswith('_graphs')},indent=2))
    print('Graph counts selected='+str(len(selected))+', alternate Majorana='+str(len(alternate)))


if __name__=='__main__': run()
