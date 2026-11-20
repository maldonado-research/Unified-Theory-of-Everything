#!/usr/bin/env python3
"""Exact literal-slot topology/derivative audit for one seed and two Y vertices.

Standard library only. Counts are labelled contractions, not amplitude weights.
No project imports or optimization-removable assertions.
"""
from collections import Counter
from fractions import Fraction as R
from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path

checks=Counter(); adverse=Counter()
def check(name,actual,expected):
    if actual!=expected:raise RuntimeError(name+': '+repr(actual)+' != '+repr(expected))
    checks[name]+=1

def reject(name,wrong,actual):
    if wrong==actual:raise RuntimeError(name+': mutation escaped')
    adverse[name]+=1

SEEDS={'S':('Lb','Lb','phi','phi'),'Sc':('L','L','phib','phib')}
YVERT={'A':('Lb','phi','nb'),'B':('L','phib','n')}
CONJ={'Lb':'L','L':'Lb','phi':'phib','phib':'phi','n':'nb','nb':'n'}
WEAK={'Lb','L','phi','phib'}
CHARGE={'Lb':R(1,2),'L':R(-1,2),'phi':R(-1,2),'phib':R(1,2),'n':R(0),'nb':R(0)}

def propagator(a,b,masses=True,wrong_charge=False):
    pair=frozenset((a,b))
    if pair==frozenset(('Lb','L')):return 'L_kinetic'
    if pair==frozenset(('phi','phib')):return 'scalar'
    if pair==frozenset(('n','nb')):return 'N_kinetic'
    if masses and a==b=='n':return 'M'
    if masses and a==b=='nb':return 'Mstar'
    if wrong_charge and a==b and a in WEAK:return 'wrong_same_charge'
    return None

def perfect_pairs(slots):
    if not slots:
        yield ();return
    first=slots[0]
    for j in range(1,len(slots)):
        for rest in perfect_pairs(slots[1:j]+slots[j+1:]):
            yield ((first,slots[j]),)+rest

def connected(edges):
    seen={0}
    while True:
        bigger=seen|{b if a in seen else a for a,b in edges if a in seen or b in seen}
        if bigger==seen:return len(seen)==3
        seen=bigger

def bridges(edges):
    return tuple(i for i in range(len(edges)) if not connected(edges[:i]+edges[i+1:]))

def enumerate_graphs(seed,ys,masses=True,wrong_charge=False):
    fields=(SEEDS[seed],)+tuple(YVERT[y] for y in ys)
    slots=tuple((v,j,f) for v,fs in enumerate(fields) for j,f in enumerate(fs))
    graphs=[]
    # Ten literal slots, six internal endpoints, three internal propagators.
    # All pairings, including potential within-vertex pairs, are considered.
    for internal in combinations(slots,6):
        for pairing in perfect_pairs(internal):
            kinds=tuple(propagator(a[2],b[2],masses,wrong_charge) for a,b in pairing)
            if None in kinds:continue
            edges=tuple((a[0],b[0]) for a,b in pairing)
            if not connected(edges):continue
            outside=tuple(s for s in slots if s not in internal)
            graphs.append({'seed':seed,'Y_vertices':ys,'internal':internal,'pairs':pairing,
                           'propagators':kinds,'edges':edges,'external':outside,
                           'bridges':bridges(edges)})
    return graphs

def sector(graph):return tuple(sorted(s[2] for s in graph['external']))
def target(seed):return tuple(sorted(SEEDS[seed]))

# A connected bridgeless three-vertex graph with three edges is a triangle.
# Include structural self-loops and parallel edges before the field test.
structural=[]
for edges in combinations_with_replacement(tuple(combinations_with_replacement(range(3),2)),3):
    degree=Counter(v for edge in edges for v in edge)
    if any(degree[v]>cap for v,cap in enumerate((4,3,3))) or not connected(edges):continue
    structural.append(edges)
proper_structure=[g for g in structural if not bridges(g)]
check('unique proper cubic-seed topology',proper_structure,[((0,1),(0,2),(1,2))])
for graph in proper_structure:
    check('proper internal degree at every vertex',tuple(sum(v in e for e in graph) for v in range(3)),(2,2,2))

# Weak Kronecker constraints can be evaluated without an assumed factor two.
# Each seed pair joins its lepton and scalar, each Yukawa joins its two
# charged slots, and every weak propagator joins its endpoint indices.
def weak_components(graph):
    vertices=(SEEDS[graph['seed']],)+tuple(YVERT[y] for y in graph['Y_vertices'])
    slots=tuple((v,j,f) for v,fs in enumerate(vertices) for j,f in enumerate(fs) if f in WEAK)
    parent={s:s for s in slots}
    def find(s):
        while parent[s]!=s:s=parent[s]
        return s
    def join(a,b):parent[find(a)]=find(b)
    for l,p in ((0,2),(1,3)):join((0,l,vertices[0][l]),(0,p,vertices[0][p]))
    for v in (1,2):join((v,0,vertices[v][0]),(v,1,vertices[v][1]))
    for a,b in graph['pairs']:
        if a[2] in WEAK:join(a,b)
    components=[]
    for root in {find(s) for s in slots}:
        components.append(tuple(s for s in slots if find(s)==root))
    return components

def weak_tensor(graph,external_values):
    external=graph['external'];weights=dict(zip(external,external_values));free=0
    for component in weak_components(graph):
        fixed={weights[s] for s in component if s in weights}
        if len(fixed)>1:return 0
        if not fixed:free+=1
    return 2**free

def derivative_data(graph):
    """Actual P0 loop coefficient for the fixed neutral pairing (0,2)(1,3).

    Orient the triangle 0->1->2->0. Incoming loop momentum at seed is
    -q at its vertex-1 line, +q at its vertex-2 line. External slots have
    zero q weight. The derivative acts on the entire second composite.
    """
    qweight={s:0 for s in graph['external']}
    for a,b in graph['pairs']:
        if a[0]==0:qweight[a]=-1 if b[0]==1 else 1
        if b[0]==0:qweight[b]=-1 if a[0]==1 else 1
    fields=SEEDS[graph['seed']]
    coefficient=qweight[(0,1,fields[1])]+qweight[(0,3,fields[3])]
    complement=qweight[(0,0,fields[0])]+qweight[(0,2,fields[2])]
    check('whole-composite loop weights are complementary',coefficient,-complement)
    seed_degree=2 if coefficient else 0
    # A direct varying external-momentum polynomial checks that this is
    # cancellation of the full composite momentum, not discarded derivatives.
    outside=[s for s in graph['external'] if s[0]==0]
    external={outside[0]:(2,1),outside[1]:(-1,3)}
    inside=[s for s in graph['internal'] if s[0]==0]
    total=tuple(sum(p[mu] for p in external.values()) for mu in range(2))
    const=dict(external)
    for slot in inside:
        const[slot]=(0,0) if qweight[slot]==-1 else tuple(-v for v in total)
    def dot(a,b):return a[0]*b[0]-a[1]*b[1]
    q=(3,1)
    def squared_polynomial(pair):
        slots=[(0,j,fields[j]) for j in pair]
        c=sum(qweight[s] for s in slots)
        k=tuple(sum(const[s][mu] for s in slots) for mu in range(2))
        return (c*c*dot(q,q),2*c*dot(q,k),dot(k,k))
    actual=squared_polynomial((1,3));other=squared_polynomial((0,2))
    check('complementary full seed momentum squares agree',actual,other)
    check('whole-composite quadratic loop polynomial',actual[0],coefficient*coefficient*8)
    if not coefficient:check('same-composite linear loop polynomial also absent',actual[1],0)
    light=sum(k=='L_kinetic' for k in graph['propagators'])
    omega=4+light+seed_degree-2*len(graph['pairs'])
    return {'seed_q_coefficient':coefficient,'seed_degree':seed_degree,
            'light_numerators':light,'superficial_UV_degree':omega}

summary=[];proper=[];by_seed={};all_graphs=[]
for seed in SEEDS:
    counts={}
    for ys in product('AB',repeat=2):
        graphs=enumerate_graphs(seed,ys);all_graphs.extend(graphs)
        accepted=[g for g in graphs if not g['bridges']]
        reducible=[g for g in graphs if g['bridges']]
        counts[''.join(ys)]={'connected':len(graphs),'proper':len(accepted),'reducible':len(reducible)}
        expected={'S':{'AA':(0,0),'AB':(8,0),'BA':(8,0),'BB':(36,12)},
                  'Sc':{'AA':(36,12),'AB':(8,0),'BA':(8,0),'BB':(0,0)}}[seed][''.join(ys)]
        check('literal connected/proper slot counts '+seed+''.join(ys),(len(graphs),len(accepted)),expected)
        for g in graphs:
            check('one-loop Euler identity',len(g['pairs'])-3+1,1)
            check('four external fields',len(g['external']),4)
            check('hypercharge conservation',sum(CHARGE[s[2]] for s in g['external']),R(0))
            if g['bridges']:
                continue
            proper.append(g)
            check('all proper graphs have no sterile external slot',any(s[2] in ('n','nb') for s in g['external']),False)
            check('proper graph has exactly one anomalous sterile propagator',sum(k in ('M','Mstar') for k in g['propagators']),1)
            check('proper graph has no kinetic sterile propagator','N_kinetic' in g['propagators'],False)
            check('proper graph contains exactly one edge of every triangle pair',set(tuple(sorted(e)) for e in g['edges']),{(0,1),(0,2),(1,2)})
            check('no free weak trace',sum(not any(s in g['external'] for s in c) for c in weak_components(g)),0)
            for values in product(range(2),repeat=4):
                check('fixed-component weak tensor is zero or one',weak_tensor(g,values) in (0,1),True)
            total_weak=sum(weak_tensor(g,values) for values in product(range(2),repeat=4))
            check('two external weak Kronecker pairs',total_weak,4)
            reject('extra weak-doublet trace factor',2*total_weak,total_weak)
            details=derivative_data(g)
            summary.append({'seed':seed,'Y_vertices':''.join(ys),'external_field_types':list(sector(g)),
                            'propagators':list(g['propagators']),**details})
    by_seed[seed]=counts
    selected=[g for g in all_graphs if g['seed']==seed and sector(g)==target(seed)]
    check('original derivative-seed target has no proper graph '+seed,sum(not g['bridges'] for g in selected),0)
    reject('connected target graphs incorrectly called proper '+seed,len(selected),0)
    # A mutation allowing absent same-charge scalar AND light propagators
    # manufactures a target triangle. The correct zero relies on both.
    wrong=sum(not g['bridges'] and sector(g)==target(seed)
              for ys in product('AB',repeat=2) for g in enumerate_graphs(seed,ys,wrong_charge=True))
    reject('same-charge charged propagators invented '+seed,wrong,0)

for seed in SEEDS:
    entries=[r for r in summary if r['seed']==seed]
    sectors=Counter(tuple(r['external_field_types']) for r in entries)
    check('four-scalar slot count '+seed,sectors[('phi','phi','phib','phib')],2)
    check('mixed LNC slot count '+seed,sectors[('L','Lb','phi','phib')],8)
    check('four-lepton slot count '+seed,sectors[('L','L','Lb','Lb')],2)
    expected=Counter({('four_scalar',2):2,('mixed',1):4,('mixed',-1):4,('four_lepton',0):2})
    def label(row):
        fs=row['external_field_types']
        return 'four_scalar' if fs==['phi','phi','phib','phib'] else 'four_lepton' if fs==['L','L','Lb','Lb'] else 'mixed'
    found=Counter((label(r),r['superficial_UV_degree']) for r in entries)
    check('whole-composite derivative degree partition '+seed,found,expected)
    reject('every mixed seed pairing falsely assigned q squared '+seed,
           sum(r['superficial_UV_degree']==1 for r in entries if label(r)=='mixed'),8)
    correct=Counter(r['superficial_UV_degree'] for r in entries)
    bad=Counter()
    for g in proper:
        if g['seed']!=seed:continue
        leg=(0,1,SEEDS[seed][1]);internal=leg in g['internal']
        bad[4+sum(k=='L_kinetic' for k in g['propagators'])+(2 if internal else 0)-6]+=1
    reject('P0 acts on lepton instead of composite '+seed,bad,correct)
    massless=sum(not g['bridges'] for ys in product('AB',repeat=2) for g in enumerate_graphs(seed,ys,masses=False))
    check('massless sterile reference has no proper triangle '+seed,massless,0)
    reject('sterile mass insertion discarded but graph retained '+seed,len(entries),massless)
    # Each target orientation has two lepton/scalar leg dressings and four
    # older CY three-point bubbles attached across the sterile bridge.
    for ys in ('AB','BA'):
        selected=[g for g in all_graphs if g['seed']==seed and ''.join(g['Y_vertices'])==ys and sector(g)==target(seed)]
        bridge_types=Counter(g['propagators'][g['bridges'][0]] for g in selected)
        check('target reducible bridge split '+seed+ys,bridge_types,
              Counter({'L_kinetic':2,'scalar':2,'N_kinetic':4}))
        reject('sterile source bridge omitted '+seed+ys,bridge_types['N_kinetic'],0)

# Independent exact chiral numerator test in Weyl blocks: P_R slash q P_R=0
# while the opposite projectors retain slash q. Gaussian rational entries.
class G:
    def __init__(self,a=0,b=0):
        if isinstance(a,G):self.a,self.b=a.a,a.b
        else:self.a,self.b=R(a),R(b)
    def __add__(self,b):b=G(b);return G(self.a+b.a,self.b+b.b)
    __radd__=__add__
    def __neg__(self):return G(-self.a,-self.b)
    def __sub__(self,b):return self+-G(b)
    def __mul__(self,b):b=G(b);return G(self.a*b.a-self.b*b.b,self.a*b.b+self.b*b.a)
    __rmul__=__mul__
    def __eq__(self,b):b=G(b);return self.a==b.a and self.b==b.b
    def conj(self):return G(self.a,-self.b)
    def __repr__(self):return 'G('+str(self.a)+','+str(self.b)+')'
def mat(n,m,f):return tuple(tuple(G(f(i,j)) for j in range(m)) for i in range(n))
def eye(n):return mat(n,n,lambda i,j:i==j)
def zero(n):return mat(n,n,lambda i,j:0)
def mul(a,b):return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(len(b))),G()) for j in range(len(b[0]))) for i in range(len(a)))
def scale(c,a):return tuple(tuple(c*v for v in row) for row in a)
def add(a,b):return tuple(tuple(x+y for x,y in zip(ar,br)) for ar,br in zip(a,b))
I=G(0,1);z2=zero(2);e2=eye(2)
pauli=(((G(),G(1)),(G(1),G())),((G(),-I),(I,G())),((G(1),G()),(G(),G(-1))))
def block(a,b,c,d):return tuple(tuple((a if i<2 and j<2 else b if i<2 else c if j<2 else d)[i%2][j%2] for j in range(4)) for i in range(4))
gamma=(block(z2,e2,e2,z2),)+tuple(block(z2,p,scale(-1,p),z2) for p in pauli)
PL=block(e2,z2,z2,z2);PR=block(z2,z2,z2,e2)
for mu in range(4):
    check('same-chirality sterile slash numerator absent',mul(mul(PR,gamma[mu]),PR),zero(4))
    check('opposite-chirality sterile slash numerator present',mul(mul(PR,gamma[mu]),PL),mul(PR,gamma[mu]))
C=scale(I,mul(gamma[2],gamma[0]));transpose=lambda a:tuple(zip(*a))
for mu in range(4):
    check('actual anomalous charge-conjugate slash vanishes',mul(mul(mul(PR,gamma[mu]),transpose(C)),transpose(PR)),zero(4))
mass_charge=mul(mul(PR,transpose(C)),transpose(PR))
reject('anomalous charge-conjugate mass numerator omitted',mass_charge,zero(4))
check('same-chirality sterile mass numerator remains',mul(PR,PR),PR)
reject('Majorana mass numerator changed to second slash',mul(PR,PR),zero(4))
M=((G(2,1),G(1,-2)),(G(1,-2),G(4,3)))
Mcal=tuple(tuple(v.conj() for v in row) for row in M)
dag=lambda a:tuple(tuple(a[j][i].conj() for j in range(len(a))) for i in range(len(a[0])))
check('left-Weyl holomorphic convention gives unbarred anomalous M',dag(Mcal),M)
reject('complex anomalous mass replaced by its conjugate',Mcal,M)

report={'status':'PASS','date':'23 October 2026 Pacific','input_commit':'f1f9a892e15545c70e4fa05827ea7fefc5a438eb',
        'calculation':'Literal one-loop 1PI topology and whole-composite derivative degree for one original seed and two pure-neutrino-Yukawa vertices',
        'counts_are':'labelled literal-slot contractions; not graph symmetry factors or amplitude multiplicities',
        'positive_checks':sum(checks.values()),'checks_by_category':dict(checks),
        'negative_checks':sum(adverse.values()),'mutations_by_category':dict(adverse),
        'connected_structural_candidates':len(structural),'proper_structural_candidates':len(proper_structure),
        'literal_graph_counts':by_seed,'proper_oriented_slot_graphs':summary,
        'restricted_target_proper_UV_pole':'zero: no allowed proper topology in this C7 Y^2 monomial',
        'scope':['unbroken zero scalar background, generic off-shell UV extraction, retained constant symmetric Majorana mass',
                 'original fixed unreduced Green/action coordinates; both seed orientations included',
                 'mixed same-composite derivative placement is UV finite; opposite-composite degree is a bound, not a pole coefficient',
                 'no integration or local coefficient claimed for the allowed four-scalar, mixed, or four-lepton descendants here',
                 'older CY source bubbles and their sterile tree bridges survive outside the proper target vertex',
                 'not complete physical C7 running, operator closure, finite matching or observable prediction'],
        'repository_imports':False,'assertions_required':False}
Path(__file__).with_name('SEED_TOPOLOGY_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
