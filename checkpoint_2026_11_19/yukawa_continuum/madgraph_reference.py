#!/usr/bin/env python3
"""Generate the chiral Yukawa UFO and evaluate the polarized three-particle cut."""
import argparse
import json
import math
from pathlib import Path
import shutil
import sys
import tempfile

from numpy.polynomial.legendre import leggauss

COMMIT = "7eba004757af0646e8e1db9d32e3d86476d9c9d5"
S = 8.0
ROOT_S = math.sqrt(S)
BEAMS = [[ROOT_S/2, 0., 0., ROOT_S/2], [ROOT_S/2, 0., 0., -ROOT_S/2]]
POINT = BEAMS + [[0.5355500529725441, -0.15655758951014173, 0.512155816530154, 0.],
                 [1.1748430830237513, -0.34344241048985824, -0.512155816530154, 0.],
                 [1.118033988749895, 0.5, 0., 0.]]


def write_ufo(path, upstream):
    path.mkdir()
    for name in ("object_library.py", "function_library.py", "coupling_orders.py"):
        shutil.copyfile(upstream/"models"/"sm"/name, path/name)

    def module(name, lines):
        (path/(name+".py")).write_text("\n".join(lines)+"\n")

    parameters = ["from object_library import all_parameters, Parameter",
                  "ZERO=Parameter(name='ZERO',nature='internal',type='real',value='0.0',texname='0')"]
    masses = [("MN1",9900011,1), ("MN2",9900012,4), ("ME1",9900013,10),
              ("ME2",9900014,20), ("MP",9900001,1)]
    for name, pdg, mass in masses:
        parameters.append(f"{name}=Parameter(name='{name}',nature='external',type='real',"
                          f"value={float(mass)},texname='{name}',lhablock='MASS',lhacode=[{pdg}])")
    module("parameters", parameters)
    particles = ["from object_library import all_particles, Particle", "import parameters as Q"]
    fields = [("e","ep","e-","e+",2,"ZERO",11,-1), ("mu","mup","mu-","mu+",2,"ZERO",13,-1),
              ("hp","hm","h+","h-",1,"MP",9900001,1)]
    for i, (mass, pdg, _) in enumerate(masses[:4], 1):
        fields.append((f"n{i}", None, f"n{i}", f"n{i}", 2, mass, pdg, 0))
    for var, conjugate, name, anti, spin, mass, pdg, charge in fields:
        particles.append(f"{var}=Particle({pdg},'{name}','{anti}',{spin},1,Q.{mass},"
                         f"Q.ZERO,'{name}','{anti}',{charge})")
        if conjugate:
            particles.append(f"{conjugate}={var}.anti()")
    module("particles", particles)
    module("lorentz", ["from object_library import all_lorentz, Lorentz",
                       "LM=Lorentz(name='LM',spins=[2,2,1],structure='ProjM(2,1)')",
                       "LP=Lorentz(name='LP',spins=[2,2,1],structure='ProjP(2,1)')"])
    couplings = ["from object_library import all_couplings, Coupling"]
    vertices = ["from object_library import all_vertices, Vertex", "import particles as P",
                "import couplings as C", "import lorentz as L"]
    for i, (lep, yukawa, conjugate) in enumerate([
            ("mu","1","1"), ("mu","0.5*complex(0,1)","-0.5*complex(0,1)"),
            ("e","1","1"), ("e","cmath.sqrt(2)*complex(0,1)","-cmath.sqrt(2)*complex(0,1)")], 1):
        pair = [(yukawa,f"P.n{i},P.{lep},P.hp","LP"), (conjugate,f"P.{lep}p,P.n{i},P.hm","LM")]
        for j, (coefficient, fields, lorentz) in enumerate(pair):
            label = 2*i+j
            couplings.append(f"G{label}=Coupling(name='G{label}',value='complex(0,1)*({coefficient})',"
                             "order={'QED':1})")
            vertices.append(f"V{label}=Vertex(name='V{label}',particles=[{fields}],color=['1'],"
                            f"lorentz=[L.{lorentz}],couplings={{(0,0):C.G{label}}})")
    module("vertices", vertices)
    module("couplings", couplings)
    objects = [("particles","particles"), ("vertices","vertices"), ("couplings","couplings"),
               ("lorentz","lorentz"), ("parameters","parameters"), ("orders","coupling_orders"),
               ("functions","function_library")]
    module("__init__", [f"import {mod}" for _, mod in objects]+
           [f"all_{name}={mod}.all_{name}" for name,mod in objects]+["gauge=[0,1]"])


def generate(upstream, model_path):
    sys.path.insert(0, str(upstream))
    from madgraph.interface.master_interface import MasterCmd
    from madgraph.core.helas_objects import HelasMatrixElement
    from madgraph.various import process_checks

    interface = MasterCmd()
    interface.exec_cmd("import model "+str(model_path), printcmd=False)
    evaluator = process_checks.MatrixElementEvaluator(interface._curr_model, reuse=False)
    matrices = {}
    counts = {}
    for label, flavor in [("full","mu"), ("pole","mu"), ("electron","e")]:
        interface.exec_cmd(f"generate {flavor}- {flavor}- > mu- n1 h-", printcmd=False)
        amplitude = interface._curr_amps[0]
        if label == "pole":
            # External legs 3 and 4 join the same intermediate charged scalar.
            def is_pole(diagram):
                return any({leg["number"] for leg in vertex["legs"][:2]} == {3,4}
                           and abs(vertex["legs"][-1]["id"]) == 9900001
                           for vertex in diagram["vertices"][:-1])
            diagrams = amplitude["diagrams"]
            amplitude["diagrams"] = type(diagrams)([dia for dia in diagrams if is_pole(dia)])
        counts[label] = len(amplitude["diagrams"])
        matrix = HelasMatrixElement(amplitude)
        generated = evaluator.evaluate_matrix_element(matrix, p=POINT, output="jamp")["jamp"]
        active = [tuple(hel) for hel, value in zip(matrix.get_helicity_matrix(), generated)
                  if abs(value[0])>1e-14]
        if set(active) != {(-1,-1,-1,-1,0), (-1,-1,-1,1,0)}:
            raise ArithmeticError("Unexpected external helicity")
        matrix_class = getattr(process_checks, "Matrix_"+matrix["processes"][0].shell_string())
        matrices[label] = matrix_class()
        matrices[label].amp2 = [0.] * counts[label]
    if counts != {"full":8, "pole":4, "electron":4}:
        raise ArithmeticError("Unexpected diagram counts")
    return matrices, evaluator.full_model, counts


def overlaps(matrices, model, momenta):
    amps = {}
    for label, matrix in matrices.items():
        matrix.clean()
        for h in (-1,1):
            matrix.matrix(momenta, [-1,-1,-1,h,0], model)
        amps[label] = [a[0] for a in matrix.jamp]
    def dot(left, right):
        return sum(x*y.conjugate() for x,y in zip(amps[left], amps[right]))
    return dot("full","electron"), dot("pole","electron"), dot("full","full")/dot("pole","pole")


def momenta(U, x, cosine, phi):
    # x is the pole scalar's beam cosine. Decay angles use its recoil axis in the LN rest frame.
    energy = (S+1-U)/(2*ROOT_S)
    momentum = math.sqrt((energy-1)*(energy+1))
    recoil = ROOT_S-energy
    rest = (U-1)/(2*math.sqrt(U))
    sintheta = math.sqrt(1-x*x)
    transverse = rest*math.sqrt(1-cosine*cosine)
    parallel = rest*(recoil*cosine-momentum)/math.sqrt(U)
    lepton = [rest*(recoil-momentum*cosine)/math.sqrt(U),
              transverse*math.cos(phi)*x+parallel*sintheta, transverse*math.sin(phi),
              -transverse*math.cos(phi)*sintheta+parallel*x]
    scalar = [energy, momentum*sintheta, 0., momentum*x]
    neutral = [ROOT_S-lepton[0]-scalar[0]]+[-lepton[j]-scalar[j] for j in (1,2,3)]
    return BEAMS+[lepton,neutral,scalar]


def electron_projection(U, x):
    h = (S-1-U)/2
    v = math.sqrt(h*h-U)
    def tree(q):
        return -30*q/((100+q)*(400+q))
    def average(d):
        return math.log1p(2*v/(d+h-v))/(2*v)
    return 2*(10*average(100)-40*average(400))/(tree(h-v*x)+tree(h+v*x))


def unit_rule(n):
    nodes, weights = leggauss(n)
    return list(zip(((1+nodes)/2).tolist(), (weights/2).tolist()))


def integrate(matrices, model, orders):
    nu, nx, nd, na = orders
    span = (ROOT_S-1)**2-1
    total = [0j,0j]
    cosines = [(2*x-1,w) for x,w in unit_rule(nd)]
    angles = [2*math.pi*(j+0.5)/na for j in range(na)]
    beam_angles = unit_rule(nx)
    for a, wa in unit_rule(nu):
        U = 1+span*(1-a*a)
        v = math.sqrt(((S-1-U)/2)**2-U)
        # (1/2) dU/(2 pi) dPhi_2(S;1,U) dPhi_2(U;0,1) with normalized angle averages.
        measure = 0.5*(2*span*a)*wa/(2*math.pi)*(2*v)/(8*math.pi*S)*(U-1)/(8*math.pi*U)
        for x, wx in beam_angles:
            scale = measure*wx*electron_projection(U,x)/na
            for cosine, wc in cosines:
                for phi in angles:
                    full, pole, _ = overlaps(matrices, model, momenta(U,x,cosine,phi))
                    total[0] += scale*wc*full
                    total[1] += scale*wc*pole
    return {"orders": orders, "full": total[0].real, "pole": total[1].real,
            "full_imaginary_residual": total[0].imag, "pole_imaginary_residual": total[1].imag,
            "ratio": (total[0]/total[1]).real}


def pair(z):
    return [z.real, z.imag]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mg5", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    upstream = args.mg5.resolve()
    if "version = 3.5.7" not in (upstream/"VERSION").read_text():
        raise RuntimeError("Expected MadGraph5_aMC@NLO 3.5.7")
    with tempfile.TemporaryDirectory() as temp:
        ufo = Path(temp)/"benchmark_yukawa"
        write_ufo(ufo, upstream)
        matrices, model, counts = generate(upstream, ufo)
        points = [POINT, momenta(1.8,0.36,0.64,0.31), momenta(2.5,0.74,-0.55,1.75)]
        samples = []
        for momenta_set in points:
            full, pole, squared = overlaps(matrices, model, momenta_set)
            samples.append({"p_mu_or_e_in_p_mu_out_p_N_p_phi": momenta_set,
                            "full_overlap": pair(full), "pole_overlap": pair(pole),
                            "overlap_ratio": pair(full/pole), "mu_square_ratio": squared.real})
        results = {"generator":"MadGraph5_aMC@NLO 3.5.7 / ALOHA (Python)", "commit": COMMIT,
                   "s_over_m2": S, "graphs": counts, "fixed_point": samples[0],
                   "other_points": samples[1:], "projected_J0": []}
        for levels in [(6,4,4,6),(10,6,6,10),(14,8,10,14),(18,12,12,16)]:
            res = integrate(matrices, model, levels)
            print("Projected J0", res, file=sys.stderr, flush=True)
            results["projected_J0"].append(res)
        coarse, fine = results["projected_J0"][-2:]
        for key in ("full","pole","ratio"):
            if not math.isclose(coarse[key], fine[key], rel_tol=1e-8):
                raise ArithmeticError("Angular or mass integration did not converge")
        args.output.write_text(json.dumps(results,indent=2)+"\n")


if __name__ == "__main__":
    main()
