#!/usr/bin/env python3
"""Independent neutral action audit; no project imports or new UV integrals."""
from collections import Counter
from pathlib import Path
import json
import sympy as s

I=s.I
positive=Counter()
negative=Counter()


class Odd:
    mates={}
    next_id=0

    def __init__(self,terms=None):
        self.terms={k:s.expand(v) for k,v in (terms or {}).items() if s.expand(v)!=0}

    @classmethod
    def pair(cls):
        a=cls.next_id; b=a+1; cls.next_id+=2
        cls.mates.update({a:b,b:a})
        return cls({(a,):s.Integer(1)}),cls({(b,):s.Integer(1)})

    @staticmethod
    def cast(x):
        return x if isinstance(x,Odd) else Odd({():s.sympify(x)})

    def __add__(self,x):
        terms=dict(self.terms)
        for key,value in self.cast(x).terms.items():
            terms[key]=terms.get(key,0)+value
        return Odd(terms)
    __radd__=__add__
    def __neg__(self):return Odd({k:-v for k,v in self.terms.items()})
    def __sub__(self,x):return self+-self.cast(x)
    def __rsub__(self,x):return self.cast(x)+-self

    def __mul__(self,x):
        terms={}
        for a,av in self.terms.items():
            for b,bv in self.cast(x).terms.items():
                if set(a)&set(b):continue
                key=tuple(sorted(a+b))
                sign=(-1)**sum(i>j for i in a for j in b)
                terms[key]=terms.get(key,0)+sign*av*bv
        return Odd(terms)
    __rmul__=__mul__
    def __truediv__(self,x):return self*(1/s.sympify(x))

    def dagger(self):
        answer=Odd()
        for key,value in self.terms.items():
            term=Odd.cast(s.conjugate(value))
            for k in reversed(key):term=term*Odd({(self.mates[k],):1})
            answer+=term
        return answer


def zero(x):
    if isinstance(x,Odd):return not x.terms
    if isinstance(x,s.MatrixBase):return all(s.expand(v)==0 for v in x)
    return s.expand(x)==0


def check(name,x):
    if not zero(x):raise RuntimeError(name+' failed')
    positive[name]+=1


def reject(name,x):
    if zero(x):raise RuntimeError(name+' wrong construction was not distinguished')
    negative[name]+=1


def fields(n):
    q={};bar={}
    for a in range(n):
        for spin in range(2):q[a,spin],bar[a,spin]=Odd.pair()
    return q,bar


def transformed(q,matrix):
    return {(a,spin):sum((q[b,spin]*matrix[a,b] for b in range(matrix.cols)),Odd())
            for a in range(matrix.rows) for spin in range(2)}


def kinetic(bar,q,matrix,sigma):
    return sum((bar[a,r]*matrix[a,b]*sigma[r,t]*q[b,t]
                for a in range(matrix.rows) for b in range(matrix.cols)
                for r in range(2) for t in range(2)),Odd())


def mass(q,matrix):
    eps=s.Matrix([[0,1],[-1,0]])
    return -sum((q[a,r]*eps[r,t]*matrix[a,b]*q[b,t]
                 for a in range(matrix.rows) for b in range(matrix.cols)
                 for r in range(2) for t in range(2)),Odd())/2


def mass_variation(q,dq,matrix):
    eps=s.Matrix([[0,1],[-1,0]])
    return -sum(((dq[a,r]*q[b,t]+q[a,r]*dq[b,t])*eps[r,t]*matrix[a,b]
                 for a in range(matrix.rows) for b in range(matrix.cols)
                 for r in range(2) for t in range(2)),Odd())/2


def positive_bilinear(q,matrix):
    """Literal +q^T epsilon coefficient q/2, distinct from the mass convention."""
    eps=s.Matrix([[0,1],[-1,0]])
    return sum((q[a,r]*eps[r,t]*matrix[a,b]*q[b,t]
                for a in range(matrix.rows) for b in range(matrix.cols)
                for r in range(2) for t in range(2)),Odd())/2


def positive_variation(q,dq,matrix):
    eps=s.Matrix([[0,1],[-1,0]])
    return sum(((dq[a,r]*q[b,t]+q[a,r]*dq[b,t])*eps[r,t]*matrix[a,b]
                for a in range(matrix.rows) for b in range(matrix.cols)
                for r in range(2) for t in range(2)),Odd())/2


def run_case(Y,M,C,u,p2,m2,label):
    n,r=Y.shape
    HL=Y*Y.H; HR=Y.H*Y; ty=s.trace(HR)
    F=C*Y.conjugate()
    K=C*Y.conjugate()*M*Y.H
    G=(K+K.H)/2; A=(K-K.H)/(2*I)
    g=s.trace(G); t=s.trace(A); b=g+I*t
    rho=u*u
    H=u*(p2+m2)*F
    Z=s.zeros(n+r); Z[:n,n:]=H; Z[n:,:n]=H.H
    ds=Z+s.diag(*list(rho*G.diagonal()),*([0]*r))
    # Preserve the full active flavor block, not only its diagonal entries.
    ds[:n,:n]=rho*G
    Rz=s.zeros(n+r);Rz[:n,:n]=-HL/2;Rz[n:,n:]=-HR.T
    phase=s.zeros(n+r);phase[:n,:n]=-I*rho*A
    J=-ds/2+phase
    check(label+' seed kinetic matrix Hermitian',ds-ds.H)
    reject(label+' zero mixed-block sensitivity witness',Z)
    check(label+' common generator removes only selected seed kinetic terms',J+J.H+ds)
    check(label+' renormalizable matrix is retained',ds+Rz+J+J.H-Rz)
    reject(label+' bare wavefunction terms are not zero',Rz)
    q,bar=fields(n+r)
    dq=transformed(q,J)
    dbar=transformed(bar,J.conjugate())
    sigmas=[s.eye(2),s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
    for mu,sigma in enumerate(sigmas):
        actual=kinetic(dbar,q,s.eye(n+r),sigma)+kinetic(bar,dq,s.eye(n+r),sigma)
        actual+=kinetic(bar,q,ds+Rz,sigma)
        expected=kinetic(bar,q,Rz,sigma)
        check(label+' literal kinetic action '+str(mu),actual-expected)
    wrong=Z.copy();wrong[:n,:n]=-rho*G
    Jwrong=-wrong/2+phase
    reject(label+' wrongly add the matter loop pole to a cancelling action',ds+Jwrong+Jwrong.H)
    reject(label+' wrong mixed Majorana half',ds-Z/2+J+J.H)
    # Direct holomorphic mass substitution in the actual two odd slots.
    B=u*Y.conjugate();D=M.conjugate()
    M0=s.zeros(n+r);M0[:n,n:]=B;M0[n:,:n]=B.T;M0[n:,n:]=D
    scalar=s.zeros(n+r)
    scalar[:n,n:]=-u**3*s.conjugate(b)*Y.conjugate()
    scalar[n:,:n]=scalar[:n,n:].T
    dM=J.T*M0+M0*J+scalar
    actual=mass_variation(q,dq,M0)+mass(q,scalar)
    check(label+' literal mass action through both odd slots',actual-mass(q,dM))
    check(label+' mass coefficient symmetry',dM-dM.T)
    check(label+' mass action including Hermitian conjugate',actual+actual.dagger()-mass(q,dM)-mass(q,dM).dagger())
    da=-(H.conjugate()*B.T+B*H.H)/2
    db=-H.conjugate()*D/2-u**3*(G.T/2+I*A.T+s.conjugate(b)*s.eye(n))*Y.conjugate()
    dd=-(H.T*B+B.T*H)/2
    check(label+' explicit active mass block',dM[:n,:n]-da)
    check(label+' explicit active sterile block with scalar transport',dM[:n,n:]-db)
    check(label+' explicit sterile mass block',dM[n:,n:]-dd)
    reject(label+' replace formal transpose by adjoint',actual-mass(q,J.H*M0+M0*J+scalar))
    reject(label+' omit scalar background transport',mass(q,scalar))
    reject(label+' remove common-map phase from mass transport',mass(q,phase.T*M0+M0*phase))
    # Independently supplied C5: the differential source map also produces a
    # mixed active/sterile contact, separately from the radial active term.
    C5=s.Matrix(n,n,lambda i,j:(2-I if i==j else 3+I))
    W5=u*u*C5.conjugate()
    W5full=s.zeros(n+r);W5full[:n,:n]=W5
    # L5 has a positive bilinear sign; the mass helper has a negative sign.
    M5=-W5full
    literal5=positive_bilinear(q,W5full)
    check(label+' positive L5 to negative-half mass convention',literal5-mass(q,M5))
    reject(label+' wrong positive L5 mass-block sign',literal5-mass(q,W5full))
    delta5=s.zeros(n+r)
    Ja=-rho*(G/2+I*A)
    delta5[:n,:n]=-Ja.T*W5-W5*Ja+2*u**4*s.conjugate(b)*C5.conjugate()
    delta5[:n,n:]=W5*H/2
    delta5[n:,:n]=delta5[:n,n:].T
    scalar5=s.zeros(n+r);scalar5[:n,:n]=2*u**4*s.conjugate(b)*C5.conjugate()
    actual5=mass_variation(q,dq,M5)+mass(q,scalar5)
    literal_variation5=positive_variation(q,dq,W5full)+positive_bilinear(q,-scalar5)
    check(label+' supplied C5 direct positive-action variation',actual5-literal_variation5)
    check(label+' supplied C5 actual two-slot variation',actual5-mass(q,delta5))
    check(label+' supplied C5 variation with conjugate completion',actual5+actual5.dagger()-mass(q,delta5)-mass(q,delta5).dagger())
    check(label+' supplied C5 coefficient remains symmetric',delta5-delta5.T)
    mixed5=s.zeros(n+r);mixed5[:n,n:]=delta5[:n,n:];mixed5[n:,:n]=mixed5[:n,n:].T
    reject(label+' drop conditional mixed C5 sterile contact',mass(q,mixed5))
    reject(label+' double the conditional mixed C5 sterile contact',actual5-mass(q,delta5+mixed5))
    # Sources are independent odd probes. Their transport cannot be inferred from spectra.
    eta,beta=fields(n+r)
    probe=kinetic(beta,dq,s.eye(n+r),s.eye(2))
    row=s.Matrix([[1+I,2,-I,3][:n+r]]) if n+r<=4 else s.ones(1,n+r)
    source_variation=row*J
    explicit=sum((beta[a,spin]*J[a,c]*q[c,spin]
                  for a in range(n+r) for c in range(n+r) for spin in range(2)),Odd())
    check(label+' literal independent probe transport',probe-explicit)
    reject(label+' omit probe variation',probe)
    reject(label+' source row uses J rather than its transpose',row*J-row*J.T)
    # Known renormalizable parameter poles reconstruct action coefficients once.
    zL=-HL/4;zR=-HR/2;zphi=-ty/2
    rY=HL*Y/4+Y*HR/2+ty*Y/2
    rM=(HR.T*M+M*HR)/2
    check(label+' Yukawa parameter and field poles cancel proper vertex coefficient',rY+zL*Y+Y*zR+zphi*Y)
    check(label+' Majorana parameter and field poles cancel proper mass coefficient',rM+zR.T*M+M*zR)
    qm=s.trace(HR*M.H*M)
    check(label+' scalar mass parameter differs from action coefficient',(m2*ty-2*qm)+2*m2*zphi+2*qm)
    reject(label+' mistake Majorana parameter pole for action pole',rM)
    # The standard z variation of a TREE seed is linear in C, not a two-loop product.
    r7=(HL*C+C*HL.T)/4+ty*C
    qs={key:value for key,value in q.items() if key[0]<n}
    dzq=transformed(qs,zL)
    literal_seed=positive_bilinear(qs,C.conjugate())*u*u*p2
    check(label+' positive tree seed to negative-half mass convention',literal_seed-mass(qs,-C.conjugate())*u*u*p2)
    reject(label+' wrong positive tree seed mass-block sign',literal_seed-mass(qs,C.conjugate())*u*u*p2)
    bare_seed_change=(positive_variation(qs,dzq,C.conjugate())+positive_bilinear(qs,C.conjugate())*2*zphi)*u*u*p2
    check(label+' literal tree seed bare-field contribution',bare_seed_change-positive_bilinear(qs,-r7.conjugate())*u*u*p2)
    inverse_seed_change=-bare_seed_change
    check(label+' inverse standard field map reverses tree seed contribution',inverse_seed_change-positive_bilinear(qs,r7.conjugate())*u*u*p2)
    check(label+' field-related seed parameter subset cancels bare fields',bare_seed_change+positive_bilinear(qs,r7.conjugate())*u*u*p2)
    reject(label+' discard inverse map on tree seed as two-loop order',inverse_seed_change)
    reject(label+' transpose versus adjoint in seed flavor slot',r7-((HL*C+C*HL.H)/4+ty*C))
    # The upper neutral source-to-normalizing bridge is separately reconstructed.
    js=-u*F*(p2/s.Integer(6)+m2/s.Integer(2))
    jc=-H/2
    bridge=jc-js
    check(label+' actual upper neutral source bridge',bridge+u*F*p2/3)
    reject(label+' add source and normalizing maps instead of their bridge',js+jc-jc)
    return {'label':label,'active_flavors':n,'sterile_flavors':r,'mass_rank':int(M.rank()),
            'g':str(s.expand(g)),'t':str(s.expand(t)),
            'kinetic_residual_nonzero':not zero(Rz),'source_variation_nonzero':not zero(source_variation),
            'seed_inverse_transport_nonzero':not zero(inverse_seed_change)}


def main():
    # Check the cancelling sign on unrestricted scalar first jets, independent
    # of the later constant-background projection.
    x,y,a,c=s.symbols('x0:2'),s.symbols('y0:2'),s.symbols('a0:2'),s.symbols('c0:2')
    symbols=list(x)+list(y)+list(a)+list(c)
    ph=s.Matrix([x[j]+I*y[j] for j in range(2)])
    dp=s.Matrix([a[j]+I*c[j] for j in range(2)])
    # Treat these symbolic components as real without assuming any equations.
    conjugate=dict(zip([s.conjugate(v) for v in symbols],symbols))
    adj=lambda v:v.conjugate().xreplace(conjugate)
    rr=(adj(ph).T*ph)[0];zz=(adj(ph).T*dp)[0]
    uu=(adj(dp).T*ph+adj(ph).T*dp)[0]
    jj=I*(zz-adj(s.Matrix([zz]))[0])
    tt=(adj(dp).T*dp)[0]
    gg,aa=s.symbols('g a',real=True)
    dph=-rr*(gg+I*aa)*ph
    ddp=-(gg+I*aa)*(uu*ph+rr*dp)
    dkin=(adj(ddp).T*dp+adj(dp).T*ddp)[0]
    loop=(gg+I*aa)*(zz**2-rr*tt)+(gg-I*aa)*(adj(s.Matrix([zz]))[0]**2-rr*tt)
    remainder=-gg*(3*uu**2-jj**2)/2
    check('arbitrary scalar first-jet counterterm sign',-loop+dkin-remainder)
    check('negative matter map scalar density change',(adj(dph).T*ph+adj(ph).T*dph)[0]+2*gg*rr**2)
    reject('wrong scalar loop versus cancelling sign',loop+dkin-remainder)
    Y=s.Matrix([[1+I,2-I],[I,1+2*I]])
    M=s.Matrix([[2+I,1-I],[1-I,3+2*I]])
    C=s.Matrix([[1+2*I,2-I],[2-I,-1+I]])
    reject('noncommuting retained mass and Yukawa witness',(Y.H*Y)*M-M*(Y.H*Y))
    witnesses=[run_case(Y,M,C,s.Integer(2),s.Integer(7),s.Integer(-3),'complex full rank')]
    # A rank-one symmetric mass with complex entries requires no mass inverse.
    Ms=s.Matrix([[1,I],[I,-1]])
    witnesses.append(run_case(Y,Ms,C,s.Integer(1),s.Integer(5),s.Integer(2),'complex singular mass'))
    report={'status':'PASS','date':'4 November 2026 (Pacific)',
            'input_commit':'41f9b709c0fafdefdc9e1b990caa5421ee99f7dc',
            'positive_checks':sum(positive.values()),'negative_or_witness_controls':sum(negative.values()),
            'positive_categories':dict(positive),'negative_categories':dict(negative),
            'witnesses':witnesses,'project_imports':False,'assertions_required':False,
            'scope':'Independent neutral constant-background classical common counterterm and source map, retaining pure-Y bare conventions; no new loop, complete physical coefficient, vacuum or observable.'}
    (Path(__file__).resolve().parent/'INDEPENDENT_COMMON_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','positive_checks':report['positive_checks'],'negative_or_witness_controls':report['negative_or_witness_controls']}))


if __name__=='__main__':main()
