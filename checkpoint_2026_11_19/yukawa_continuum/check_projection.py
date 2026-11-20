"""Project the finite scalar insertion and full explicit scale log onto rate observables."""
import argparse
import json
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / "radiative_stability"
ENERGIES = (6, 8, 10)


def load(path):
    return json.loads(path.read_text(), parse_float=Decimal)


def determinant(a, b):
    """Determinant with columns (1,a,b)."""
    return (a[1]-a[0])*(b[2]-b[0]) - (a[2]-a[0])*(b[1]-b[0])


def close(value, expected, relative, label):
    if abs(value-expected) > relative*abs(expected):
        raise ArithmeticError(label)


def project(ell, d, e, endpoint_ratio):
    det = determinant(d, e)
    weight = (d[1]-d[0])/(d[2]-d[0])
    contrast = ell[1]-ell[0]-weight*(ell[2]-ell[0])
    closure = determinant(d, ell)/det
    if contrast != -det*closure/(d[2]-d[0]):
        raise ArithmeticError("the two shape projections must agree up to their normalizations")
    # Solve all three normalized rates for a finite change h in the linear
    # model. The central difference tests differentiation of the nonlinear
    # invariant itself, including terms from the fitted normalization.
    def inferred_invariant(h):
        observed = tuple(1+h*v for v in ell)
        W = determinant(d, observed)/det
        V = determinant(observed, e)/det
        U = observed[0]-d[0]*V-e[0]*W
        if U <= 0:
            raise ArithmeticError("positive normalization required")
        return W/U-(V/U)**2
    h = Fraction(1, 10**7)
    central = (inferred_invariant(h)-inferred_invariant(-h))/(2*h)
    close(central, closure, Fraction(1, 10**8), "nonlinear closure derivative")
    return {"fractional_rate_shifts":ell, "endpoint_ratio_shift":endpoint_ratio*(ell[2]-ell[0]),
            "log_contrast_shift":contrast,"four_rate_invariant_shift":closure}


def serial(projected):
    return {name:([float(f"{float(t):.9g}") for t in value] if isinstance(value, tuple)
                  else float(f"{float(value):.9g}")) for name,value in projected.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    detuning=load(PRIOR/"DETUNING.json")["s_over_m2"]
    d=tuple(Fraction(detuning[str(s)]["dlog_rate_dRe_delta_at_zero"]) for s in ENERGIES)
    e=tuple(Fraction(detuning[str(s)]["c"])/Fraction(detuning[str(s)]["a"]) for s in ENERGIES)
    quartic=load(PRIOR/"RESCATTERING.json")
    ratio=Fraction(quartic["leading_rate_ratio_10_over_6"])
    reference=project(tuple(Fraction(t["fractional_rate_derivative_at_lambda_zero"])
                            for t in quartic["rows"]), d,e,ratio)
    old=load(PRIOR/"CLOSURE.json")["quartic_solution_derivatives_U_V_W_over_normalization"][2]
    close(reference["four_rate_invariant_shift"], Fraction(old), Fraction(1,10**9),"quartic control")
    simplex=load(HERE/"CONTINUUM.json")["s_over_m2"]
    external=load(HERE/"QCDLOOP_CONTINUUM.json")["rows"]
    # The external library receipt lists successive quadrature refinements.
    fine={s:max((row for row in external if row["s_over_m2"]==s),
                key=lambda row:(row["mass_order"],row["angle_order"])) for s in ENERGIES}
    logs=load(HERE/"RG_LOG.json")["rows"]
    if [row["s_over_m2"] for row in logs]!=list(ENERGIES):
        raise ArithmeticError("unexpected scale-log energy set")
    finite,scale_log={},{}
    for block,flavor in (("g","muon"),("y","electron")):
        def ell(row):
            return Fraction(row["m2_sigma_slope_"+block])/Fraction(row["m2_sigma0"])
        primary=project(tuple(ell(simplex[str(s)]) for s in ENERGIES),d,e,ratio)
        secondary=project(tuple(ell(fine[s]) for s in ENERGIES),d,e,ratio)
        for key in ("endpoint_ratio_shift","log_contrast_shift","four_rate_invariant_shift"):
            close(primary[key], secondary[key], Fraction(1,10**5),"simplex vs QCDLoop for "+key)
        finite[flavor]=serial(primary)
        log=tuple(Fraction(row["fractional_explicit_log_per_"+flavor]) for row in logs)
        scale_log[flavor]=serial(project(log,d,e,ratio))
    data={"energy_order_s_over_m2":ENERGIES, "quartic_reference":serial(reference),
          "finite_scalar_insertion_per_squared_base_yukawa":finite,
          "explicit_scale_log_coefficient_per_squared_base_yukawa":scale_log}
    text=json.dumps(data,indent=2,allow_nan=False)+"\n"
    if args.output:
        args.output.write_text(text)
    else:
        print(text,end="")


if __name__=="__main__":
    main()
