#!/usr/bin/env python3
"""Exact engineering, spurion and bare-mu grading; no loop integration."""
from collections import Counter
from dataclasses import dataclass
import json
from pathlib import Path
import sympy as s


positive=Counter()
adverse=Counter()


def check(name,value):
    if s.simplify(value)!=0:
        raise RuntimeError(name+' failed: '+str(value))
    positive[name]+=1


def reject(name,value):
    if s.simplify(value)==0:
        raise RuntimeError(name+' incorrect control escaped')
    adverse[name]+=1


# Four-dimensional mass dimension, explicit mu^(weight*epsilon), lepton charge.
PARAMETERS={
    'C7':(-3,2,2), 'C7c':(-3,2,-2), 'C5':(-1,2,2),
    'Y':(0,1,0), 'Yc':(0,1,0),
    'M':(1,0,-2), 'Md':(1,0,2),
    'm2':(2,0,0), 'lambda':(0,2,0),
}
# Type and global lepton charge. Scalars include their conjugates as distinct slots.
FIELDS={
    'phi':('s',0), 'phid':('s',0),
    'L':('f',1), 'barL':('f',-1), 'Lc':('f',-1),
    'NL':('f',-1), 'barNL':('f',1),
    'NR':('f',1), 'barNR':('f',-1),
}


@dataclass(frozen=True)
class Row:
    name:str
    fields:tuple
    derivatives:int
    parameters:tuple
    description:str

    def ledger(self):
        ns=sum(FIELDS[f][0]=='s' for f in self.fields)
        nf=sum(FIELDS[f][0]=='f' for f in self.fields)
        dimension=s.Integer(ns)+s.Rational(3,2)*nf+self.derivatives
        charge=sum(FIELDS[f][1] for f in self.fields)
        p_dim=sum(PARAMETERS[p][0] for p in self.parameters)
        p_weight=sum(PARAMETERS[p][1] for p in self.parameters)
        p_charge=sum(PARAMETERS[p][2] for p in self.parameters)
        return {
            'name':self.name,'description':self.description,
            'scalar_slots':ns,'fermion_slots':nf,'derivatives':self.derivatives,
            'operator_dimension_4d':int(dimension),
            'coefficient_dimension_4d':p_dim,
            'recipient_mu_weight':ns+nf-2,
            'input_mu_weight':p_weight,'weight_difference':p_weight-(ns+nf-2),
            'operator_lepton_charge':charge,'coefficient_lepton_charge':p_charge,
            'Yukawa_degree':sum(p in {'Y','Yc'} for p in self.parameters),
            'sterile_mass_degree':sum(p in {'M','Md'} for p in self.parameters),
            'C7_degree':sum(p in {'C7','C7c'} for p in self.parameters),
            'C5_degree':self.parameters.count('C5'),
            'mass_squared_degree':self.parameters.count('m2'),
            'quartic_degree':self.parameters.count('lambda'),
            'fixed_mass_inverse_scale_degree':3*sum(p in {'C7','C7c'} for p in self.parameters)
                                               +self.parameters.count('C5'),
            'parameter_tokens':list(self.parameters),
        }


ROWS=[
    Row('seed external-field input',('barL','Lc','phi','phi'),2,
        ('C7','Y','Yc'),'Known field contribution in the original seed coordinates.'),
    Row('three-point derivative source',('barL','NL','phi'),3,
        ('C7','Yc'),'Derivative-three part of the earlier F=C7 Y* counterterm.'),
    Row('three-point mass-dressed source',('barL','NL','phi'),1,
        ('C7','Yc','m2'),'The m_H^2 F term has a dimension-five field monomial.'),
    Row('four-lepton primitive',('barL','barL','L','L'),0,
        ('C7','Yc','M','Yc'),'Displayed orientation of the four-lepton mass pole.'),
    Row('four-scalar derivative primitive',('phi','phi','phid','phid'),2,
        ('C7','Yc','M','Yc'),'Real and imaginary scalar derivative placements.'),
    Row('mixed derivative primitive',('barL','L','phi','phid'),1,
        ('C7','Yc','M','Yc'),'Mixed two-lepton/two-scalar mass pole.'),
    Row('quartic mass-cubic primitive',('phi','phi','phid','phid'),0,
        ('C7','Yc','M','Md','M','Yc'),'Ordered M M^dagger M has net charge -2.'),
    Row('radial mass-squared potential transport',('phi','phi','phid','phid'),0,
        ('C7','Yc','M','Yc','m2'),'Potential change from the radial normalization.'),
    Row('induced neutrino Yukawa transport',('barL','NR','phi','phi','phid'),0,
        ('C7','Yc','M','Yc','Y'),'rho times the declared neutrino Yukawa interaction.'),
    Row('induced scalar sextic transport',('phi','phi','phi','phid','phid','phid'),0,
        ('C7','Yc','M','Yc','lambda'),'rho^3 term induced by the reference quartic.'),
    Row('conditional dimension-five transport',('barL','Lc','phi','phi','phi','phid'),0,
        ('C5','C7','Yc','M','Yc'),'Requires an independently supplied symmetric C5.'),
    Row('conditional C5 source-map derivative transport',('barNL','Lc','phi','phi','phid'),2,
        ('C5','C7c','Y'),'Derivative part of C5 F* P, with F=C7 Y*.'),
    Row('conditional C5 source-map mass-squared transport',('barNL','Lc','phi','phi','phid'),0,
        ('C5','C7c','Y','m2'),'Mass-squared part of C5 F* P, P=-Box+m_H^2.'),
    Row('renormalizable lepton field input',('barL','L'),1,
        ('Y','Yc'),'Pure-neutrino-Yukawa lepton field pole.'),
    Row('renormalizable scalar field input',('phi','phid'),2,
        ('Y','Yc'),'Pure-neutrino-Yukawa scalar field pole.'),
    Row('renormalizable Majorana mass input',('barNL','NR'),0,
        ('M','Y','Yc'),'One displayed mass-spurion orientation of the mass input.'),
    Row('renormalizable Yukawa parameter input',('barL','NR','phi'),0,
        ('Y','Yc','Y'),'Degree-three Yukawa parameter/field combination.'),
    Row('quartic field-factor transport',('phi','phi','phid','phid'),0,
        ('lambda','Y','Yc'),'Reference quartic varied by its known scalar field pole.'),
]


def signature(row):
    """A necessary grading test, not a complete independent-operator quotient."""
    item=row.ledger()
    return (item['scalar_slots'],item['fermion_slots'],item['derivatives'],
            item['operator_dimension_4d'],item['operator_lepton_charge'])


def main():
    ns,nf,nd,eps=s.symbols('ns nf nd eps')
    d_operator=ns*(1-eps)+nf*(s.Rational(3,2)-eps)+nd
    d_recipient=4-2*eps-d_operator
    check('general d-dimensional coefficient dimension',
          d_recipient-(4-ns-s.Rational(3,2)*nf-nd+(ns+nf-2)*eps))
    external,internal,vertices,loops=s.symbols('E I V loops')
    total_vertex_slots=external+2*internal
    weight_difference=(total_vertex_slots-2*vertices)-(external-2)
    check('connected graph canonical-weight identity',
          weight_difference.subs(internal,loops+vertices-1)-2*loops)
    win,wout,r,beta=s.symbols('win wout r beta')
    # At fixed bare coefficient the finite one-loop condition is
    # beta + w_out*r - w_in*r = 0, after canonical -epsilon weights.
    beta_solution=s.solve(beta+wout*r-win*r,beta)[0]
    check('signed parameter-pole weight conversion',beta_solution-(win-wout)*r)

    # Literal reference actions fix the spurion signs independently of the rows.
    reference=[
        ('Yukawa',('barL','NR','phi'),0,('Y',),4,1),
        ('Majorana mass',('barNL','NR'),0,('M',),3,0),
        ('dimension-five input',('barL','Lc','phi','phi'),0,('C5',),5,2),
        ('dimension-seven seed',('barL','Lc','phi','phi'),2,('C7',),7,2),
        ('scalar quartic',('phi','phi','phid','phid'),0,('lambda',),4,2),
    ]
    for name,fields,derivatives,parameters,dimension,mu_weight in reference:
        data=Row(name,fields,derivatives,parameters,'reference action').ledger()
        check('reference dimension '+name,data['operator_dimension_4d']-dimension)
        check('reference recipient mu '+name,data['recipient_mu_weight']-mu_weight)
        check('reference coefficient dimension '+name,
              data['coefficient_dimension_4d']+dimension-4)
        check('reference spurion invariance '+name,
              data['operator_lepton_charge']+data['coefficient_lepton_charge'])
        check('reference bare mu consistency '+name,
              data['input_mu_weight']-data['recipient_mu_weight'])

    ledger=[]
    for row in ROWS:
        data=row.ledger();ledger.append(data)
        check('engineering dimension '+row.name,
              data['operator_dimension_4d']+data['coefficient_dimension_4d']-4)
        check('spurion invariance '+row.name,
              data['operator_lepton_charge']+data['coefficient_lepton_charge'])
        check('one-loop weight difference '+row.name,data['weight_difference']-2)
    by_name={row.name:row for row in ROWS}
    # Vertex valences are topological consistency checks, not existence proofs.
    for name,valences,external_slots in [
        ('seed plus one Yukawa',(4,3),3),
        ('seed plus two Yukawas',(4,3,3),4),
        ('one-loop scalar two-point',(3,3),2),
        ('seed, two Yukawas and quartic',(4,3,3,4),6),
    ]:
        internal_lines=(sum(valences)-external_slots)//2
        number_loops=internal_lines-len(valences)+1
        check('valence has integer internal lines '+name,
              sum(valences)-external_slots-2*internal_lines)
        check('valence one-loop consistency '+name,number_loops-1)
        check('valence canonical-weight conversion '+name,
              sum(n-2 for n in valences)-(external_slots-2)-2)

    reject('use weight two for the three-point recipient',
           by_name['three-point derivative source'].ledger()['input_mu_weight']-2-2)
    reject('use weight two for the scalar sextic recipient',
           by_name['induced scalar sextic transport'].ledger()['input_mu_weight']-2-2)
    reject('use weight two for the induced Yukawa recipient',
           by_name['induced neutrino Yukawa transport'].ledger()['input_mu_weight']-2-2)
    reject('use weight two for the conditional dimension-seven recipient',
           by_name['conditional dimension-five transport'].ledger()['input_mu_weight']-2-2)
    wrong_mass=Row('wrong literal M cubed',('phi','phi','phid','phid'),0,
                   ('C7','Yc','M','M','M','Yc'),'incorrect complex mass word').ledger()
    reject('replace M Mdag M by literal M cubed in the spurion grading',
           wrong_mass['coefficient_lepton_charge']+wrong_mass['operator_lepton_charge'])
    missing_mass=Row('missing retained mass',('barL','L','phi','phid'),1,
                    ('C7','Yc','Yc'),'incorrect omission').ledger()
    reject('omit the retained mass from the lower-dimensional mixed pole',
           missing_mass['operator_dimension_4d']+missing_mass['coefficient_dimension_4d']-4)
    reject('call the mass-dressed source field monomial dimension seven',
           by_name['three-point mass-dressed source'].ledger()['operator_dimension_4d']-7)
    reject('give derivatives a scalar epsilon dimension shift',
           by_name['three-point derivative source'].derivatives)
    reject('identify scalar and four-lepton terms from dimension and weight alone',
           int(signature(by_name['four-lepton primitive'])!=signature(by_name['four-scalar derivative primitive'])))
    reject('identify the seed and conditional C5 transport from dimension alone',
           int(signature(by_name['seed external-field input'])!=signature(by_name['conditional dimension-five transport'])))
    kappa=s.symbols('kappa')
    loop_order=s.degree(kappa*kappa,kappa)
    reject('retain a product of two one-loop directions at first loop order',loop_order-1)
    # Separate inverse-scale and loop degrees; neither is a substitute for the other.
    check('seed-linear plus seed-linear product is seed quadratic',
          by_name['three-point derivative source'].ledger()['C7_degree']
          +by_name['mixed derivative primitive'].ledger()['C7_degree']-2)
    check('conditional independent input is linear in each coefficient',
          by_name['conditional dimension-five transport'].ledger()['C7_degree']-1)
    check('conditional input has one separate C5 factor',
          by_name['conditional dimension-five transport'].ledger()['C5_degree']-1)
    check('canonical mass degree differs from mass-spurion charge',
          by_name['quartic mass-cubic primitive'].ledger()['sterile_mass_degree']-3)
    check('conditional matter C5 term has inverse-scale degree four',
          by_name['conditional dimension-five transport'].ledger()['fixed_mass_inverse_scale_degree']-4)
    check('conditional source-map C5 derivative term has inverse-scale degree four',
          by_name['conditional C5 source-map derivative transport'].ledger()['fixed_mass_inverse_scale_degree']-4)
    check('conditional source-map mass-squared term has inverse-scale degree four',
          by_name['conditional C5 source-map mass-squared transport'].ledger()['fixed_mass_inverse_scale_degree']-4)
    reject('retain a conditional C5C7 term in a strict inverse-scale-degree-three truncation',
           by_name['conditional dimension-five transport'].ledger()['fixed_mass_inverse_scale_degree']-3)

    report={
        'status':'PASS','date':'4 November 2026 (Pacific)',
        'input_commit':'41f9b709c0fafdefdc9e1b990caa5421ee99f7dc',
        'scope':'Exact engineering/spurion/bare-mu grading of specified prior poles and action transports; no new loop, graph-existence proof, operator-basis quotient or total physical beta function.',
        'positive_checks':sum(positive.values()),'rejected_adverse_controls':sum(adverse.values()),
        'positive_categories':dict(positive),'adverse_categories':dict(adverse),
        'general_identities':{
            'recipient_mu_weight':'n_scalar+n_fermion-2 in the declared zero-gauge sector',
            'connected_graph_weight_difference':'sum_v(N_v-2)-(E-2)=2L',
            'signed_parameter_residue_conversion':'16 pi^2 beta_c|(signed residue r)=(w_in-w_out)r',
            'fixed_mass_inverse_scale_degree':'3*n_C7+n_C5 under C7~Lambda^-3, C5~Lambda^-1 and fixed retained masses',
        },
        'rows':ledger,'project_imports':False,'assertions_required':False,
    }
    (Path(__file__).resolve().parent/'SECTOR_GRADING_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','rows':len(ledger),'positive_checks':report['positive_checks'],
                      'rejected_adverse_controls':report['rejected_adverse_controls']}))


if __name__=='__main__':
    main()
