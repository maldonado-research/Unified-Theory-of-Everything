#!/usr/bin/env python3
"""Independent SU(2) index and one-leg UV audit; standard library only.

No repository imports. Triangles are integrated by their asymptotic angular
UV functional, independently of the displayed R_L/R_H formulas. The present
checker does not perform an EOM reduction, a two-W match, or physical RGE.
"""
from fractions import Fraction as Q
from itertools import product
from math import prod
import json

ETA=(1,-1,-1,-1)
ZERO=(0,0,0,0)
COUNTS={}

def check(a,b,label):
    if a!=b: raise RuntimeError(label+": "+repr((a,b)))
    COUNTS[label.split(':')[0]]=COUNTS.get(label.split(':')[0],0)+1

def reject(a,b,label):
    if a==b: raise RuntimeError('mutation survived: '+label)
    COUNTS['mutation_rejections']=COUNTS.get('mutation_rejections',0)+1

def clean(a): return {k:v for k,v in a.items() if v}
def ca(*aa):
    z={}
    for a in aa:
        for k,v in a.items(): z[k]=z.get(k,Q(0))+v
    return clean(z)
def cs(c,a): return clean({k:c*v for k,v in a.items()})
def cm(a,b):
    z={}
    for i,x in a.items():
        for j,y in b.items():
            swaps=sum((j&((1<<r)-1)).bit_count() for r in range(4) if i&(1<<r))
            sign=(-1)**swaps*prod(ETA[r] for r in range(4) if (i&j)&(1<<r))
            k=i^j
            z[k]=z.get(k,Q(0))+sign*x*y
    return clean(z)
def chain(*xs):
    z={0:Q(1)}
    for x in xs: z=cm(z,x)
    return z
def slash(v): return clean({1<<i:Q(ETA[i])*v[i] for i in range(4)})
def dot(a,b): return sum(Q(ETA[i])*a[i]*b[i] for i in range(4))
def va(*vv): return tuple(sum(v[i] for v in vv) for i in range(4))
def vs(c,v): return tuple(c*x for x in v)

def pa(*pp):
    z={}
    for p in pp:
        for k,v in p.items(): z[k]=ca(z.get(k,{}),v)
    return clean(z)
def ps(c,p): return clean({k:cs(c,v) for k,v in p.items()})
def pm(a,b):
    z={}
    for i,x in a.items():
        for j,y in b.items():
            k=tuple(i[r]+j[r] for r in range(4))
            z[k]=ca(z.get(k,{}),cm(x,y))
    return clean(z)
def pc(x): return {} if not x else {ZERO:x}
def pn(x): return pc({0:Q(x)})
def qcomponent(r):
    e=tuple(int(i==r) for i in range(4))
    return {e:{0:Q(1)}}
QV=tuple(qcomponent(i) for i in range(4))
QSL=pa(*(ps(ETA[i],pm(QV[i],pc({1<<i:Q(1)}))) for i in range(4)))
QSQ=pa(*(ps(ETA[i],pm(QV[i],QV[i])) for i in range(4)))
def qdot(v): return pa(*(ps(ETA[i]*v[i],QV[i]) for i in range(4)))
def qsqshift(v): return pa(QSQ,ps(-2,qdot(v)),pn(dot(v,v)))
def homogeneous(p,n): return {k:v for k,v in p.items() if sum(k)==n}
def angle(p):
    z={}
    for powers,c in p.items():
        if any(x%2 for x in powers): continue
        n=sum(powers)//2
        factor=Q(prod(prod(range(1,x,2))*ETA[i]**(x//2) for i,x in enumerate(powers)),prod(4+2*j for j in range(n)))
        z=ca(z,cs(factor,c))
    return z
def triangle(numerator,shifts,masses2):
    """Logarithmic UV coefficient from large-q expansion at generic momenta.

    q² terms need no expansion; q³ needs 2q.sum(a); q⁴ needs
    -sum(a²-m²)+4 sum_{i<=j}(q.a_i)(q.a_j). The angular moments are
    g/4, sym(gg)/24 and sym(ggg)/192. This defines the UV functional,
    not a prescription to discard scaleless original massless integrals.
    """
    n2,n3,n4=(homogeneous(numerator,n) for n in (2,3,4))
    shift_sum=va(*shifts)
    constant=sum(dot(a,a)-m for a,m in zip(shifts,masses2))
    pairs=pa(*(pm(qdot(a),qdot(b)) for i,a in enumerate(shifts) for b in shifts[i:]))
    return ca(angle(n2),cs(2,angle(pm(n3,qdot(shift_sum)))),
              cs(-constant,angle(n4)),cs(4,angle(pm(n4,pairs))))

def integrate(p,h,k,e,m2,direct=False):
    t=va(p,vs(-1,k)); v=va(h,k); u=va(vs(2,p),vs(-1,k))
    seed_l=pn(2*dot(t,t)) if direct else qsqshift(v)
    seed_h=pn(2*dot(t,t)) if direct else qsqshift(h)
    fermion=pm(pm(pa(QSL,pc(cs(-1,slash(k)))),pc(slash(e))),QSL)
    scalar=pm(QSL,pa(pn(dot(u,e)),ps(-2,qdot(e))))
    return (triangle(pm(seed_l,fermion),(ZERO,k,p),(0,0,m2)),
            triangle(pm(seed_h,scalar),(ZERO,p,t),(0,m2,m2)))

def green(p,h,k,e,m2):
    r=va(p,vs(-1,h),vs(-1,k))
    return ca(cs(-dot(va(vs(2,r),k),e)/3+dot(va(vs(2,h),k),e)/6,slash(p)),
              cs(dot(p,p)/6-m2/2,slash(e)))
def transverse(p,h,k,e):
    r=va(p,vs(-1,h),vs(-1,k))
    comm=ca(chain(slash(k),slash(e)),cs(-1,chain(slash(e),slash(k))))
    return ca(cs(Q(1,4),chain(comm,slash(r))),
              cs(Q(1,3),ca(cs(dot(k,e),slash(k)),cs(-dot(k,k),slash(e)))))
def bubble(p,h,m2):
    return ca(cs(dot(p,p)+m2-Q(2,3)*dot(p,h)+dot(h,h)/2,slash(p)),
              cs(dot(p,p)/6-m2/2,slash(h)))

# Exact complex 2x2 weak matrices represented by rational real/imag pairs.
def z(x=0,y=0): return (Q(x),Q(y))
def za(*xx): return (sum(x[0] for x in xx),sum(x[1] for x in xx))
def zm(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def zs(c,a): return (c*a[0],c*a[1])
def mat(rows): return tuple(tuple(x if isinstance(x,tuple) else z(x) for x in row) for row in rows)
def ma(*aa): return tuple(tuple(za(*(a[i][j] for a in aa)) for j in range(2)) for i in range(2))
def ms(c,a): return tuple(tuple(zs(c,x) for x in row) for row in a)
def mt(a): return tuple(tuple(a[j][i] for j in range(2)) for i in range(2))
def mm(a,b): return tuple(tuple(za(*(zm(a[i][j],b[j][k]) for j in range(2))) for k in range(2)) for i in range(2))
EPS=mat([[0,1],[-1,0]])
ID=mat([[1,0],[0,1]])
ZM=ms(0,ID)
TS=(ms(Q(1,2),mat([[0,1],[1,0]])),ms(Q(1,2),mat([[0,z(0,-1)],[z(0,1),0]])),ms(Q(1,2),mat([[1,0],[0,-1]])))

def phi_index_route(t):
    tb=ms(-1,mt(t))
    # delta_ai delta_bj is direct; delta_aj delta_bi is crossed.
    # The Yukawa contracts charge-conjugate fermion e with scalar f.
    out=[]
    for cross,fermion in product((False,True),repeat=2):
        rows=[]
        for a in range(2):
            row=[]
            for c in range(2):
                terms=[]
                for b,d,e in product(range(2),repeat=3):
                    seed=(int(a==d)*int(b==c)) if cross else (int(a==c)*int(b==d))
                    insertion=tb[b][e] if fermion else t[d][e]
                    other=int(d==e) if fermion else int(b==e)
                    terms.append(zs(seed*other,insertion))
                row.append(za(*terms))
            rows.append(tuple(row))
        out.append(tuple(rows))
    return out

def epsilon_index_route(t):
    tb=ms(-1,mt(t)); out=[]
    for cross,fermion in product((False,True),repeat=2):
        rows=[]
        for a in range(2):
            row=[]
            for c in range(2):
                terms=[]
                for b,d,e in product(range(2),repeat=3):
                    seed=zm(EPS[a][d],EPS[b][c]) if cross else zm(EPS[a][c],EPS[b][d])
                    route=zm(tb[b][e],EPS[e][d]) if fermion else zm(tb[d][e],EPS[b][e])
                    terms.append(zm(seed,route))
                row.append(za(*terms))
            rows.append(tuple(row))
        out.append(tuple(rows))
    return out

def weak_checks():
    for i,t in enumerate(TS):
        check(ma(mm(t,EPS),mm(EPS,mt(t))),ZM,'epsilon_pseudoreality:'+str(i))
        check(phi_index_route(t),[ZM,ZM,t,ms(-1,t)],'phi_group_route:'+str(i))
        te=mm(t,EPS)
        check(epsilon_index_route(t),[ZM,ZM,te,ms(-1,te)],'epsilon_group_route:'+str(i))
        reject(phi_index_route(t)[0],ms(2,t),'direct factor 2T is false')
        reject(phi_index_route(t)[3],t,'conjugate-fermion crossed sign')
    hyper=ms(-Q(1,2),ID)
    check(phi_index_route(hyper),[ms(2,hyper),ms(-2,hyper),hyper,ms(-1,hyper)],'hypercharge_control')
    q=ma(hyper,TS[2])
    check(q,mat([[0,0],[0,-1]]),'photon_projection')
    reject(hyper[0][0],q[0][0],'B-only neutral photon projection')
    # Exact rational angle is a check of the algebraic formulas in the note.
    sw,cw,ee=Q(3,5),Q(4,5),Q(3,7)
    gw,gp=ee/sw,ee/cw
    photon=ma(ms(cw*gp,hyper),ms(sw*gw,TS[2]))
    zfield=ma(ms(-sw*gp,hyper),ms(cw*gw,TS[2]))
    check(photon,ms(ee,q),'photon_projection_with_couplings')
    check(zfield[0][0],z(ee/(2*sw*cw)),'neutral_Z_projection')
    reject(zfield[0][0],z(),'neutral_Z_nonzero')
    check(photon[1][1],z(-ee),'charged_photon_control')

def case(p,h,k,e,m2,label):
    fl,sh=integrate(p,h,k,e,m2)
    loop=ca(sh,cs(-1,fl))
    expect=ca(green(p,h,k,e,m2),transverse(p,h,k,e))
    check(loop,expect,'crossed_pole:'+label)
    dl,ds=integrate(p,h,k,e,m2,True)
    t=va(p,vs(-1,k))
    expectdirect=cs(-dot(t,t),slash(e))
    check(dl,expectdirect,'direct_fermion_pole:'+label)
    check(ds,expectdirect,'direct_scalar_pole:'+label)
    check(ca(ds,cs(-1,dl)),{},'direct_pole_cancellation:'+label)
    if any(k):
        flk,shk=integrate(p,h,k,k,m2)
        check(ca(shk,cs(-1,flk)),ca(bubble(p,va(h,k),m2),cs(-1,bubble(p,h,m2))),
              'ward_identity:'+label)
    return loop

def main():
    weak_checks()
    bases=[tuple(Q(int(i==j)) for i in range(4)) for j in range(4)]
    # Degree-two momentum polynomial: diagonal and pairwise sum samples
    # determine all coefficients; polarization enters linearly. Mass² is
    # independent and degree one. Each sample is of the local UV polynomial.
    for eidx,e in enumerate(bases):
        for i in range(12):
            for j in range(i,12):
                vv=[Q(int(a==i)+int(a==j and j!=i)) for a in range(12)]
                case(tuple(vv[:4]),tuple(vv[4:8]),tuple(vv[8:]),e,Q(0),f'{eidx},{i},{j}')
        case(ZERO,ZERO,ZERO,e,Q(1),f'mass,{eidx}')
    p=tuple(map(Q,(3,-2,1,4)));h=tuple(map(Q,(1,2,-3,1)));k=tuple(map(Q,(2,1,-1,-2)));e=tuple(map(Q,(1,-1,2,3)))
    loop=case(p,h,k,e,Q(-5,3),'generic')
    reject(loop,green(p,h,k,e,Q(-5,3)),'omitted transverse polynomial')
    fl,sh=integrate(p,h,k,e,Q(-5,3))
    reject(loop,ca(sh,fl),'equal effective crossed charges')
    report={'status':'PASS','counts':COUNTS,'total_checks':sum(COUNTS.values()),
            'result':'SU(2) direct trace vanishes; crossed +/-T; direct abelian UV poles cancel; R_W=G_J+D; neutral photon pole cancels',
            'limits':'one external W only; no physical anomalous dimension; no two-W or mixed BW matching; no EOM reduction'}
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=='__main__': main()
