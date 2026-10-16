#!/usr/bin/env python3
"""Verify the public payload and run checks on a disposable copy."""
from pathlib import Path, PurePosixPath
import argparse, hashlib, json, platform, shutil, subprocess, sys, tempfile
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parent


def verify_manifest():
    manifest=json.loads((ROOT/'PUBLIC_MANIFEST.json').read_text())
    seen=set()
    for row in manifest['files']:
        rel=row['path'];path=PurePosixPath(rel)
        if path.is_absolute() or '..' in path.parts or rel in seen:
            raise RuntimeError('Unsafe or repeated manifest path: '+rel)
        seen.add(rel);target=ROOT/rel
        if target.is_symlink() or not target.is_file():
            raise RuntimeError('Missing or nonregular payload: '+rel)
        data=target.read_bytes()
        if len(data)!=row['bytes'] or hashlib.sha256(data).hexdigest()!=row['sha256']:
            raise RuntimeError('Payload hash mismatch: '+rel)
    actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
            if p.is_file() and '.git' not in p.relative_to(ROOT).parts
            and '.venv' not in p.relative_to(ROOT).parts
            and '__pycache__' not in p.relative_to(ROOT).parts
            and p.name not in {'.DS_Store'}
            and p.relative_to(ROOT).as_posix() not in manifest['excluded_self_and_receipts']}
    if actual!=seen:
        raise RuntimeError('Unexpected/missing payload: '+str(sorted(actual.symmetric_difference(seen))))
    return len(seen)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full',action='store_true',help='Include the numerical pulse benchmark and separate Radau comparison')
    parser.add_argument('--output',type=Path,help='Optional replay receipt destination; must be outside the repository')
    args=parser.parse_args()
    if sys.flags.optimize:
        raise SystemExit('Run verify.py without -O; checkpoint review assertions must remain active.')
    if args.output and (args.output.resolve()==ROOT or ROOT in args.output.resolve().parents):
        raise SystemExit('The replay receipt destination must be outside the repository.')
    count=verify_manifest();print(f'Public payload: PASS ({count} files)',flush=True)
    results=[]
    with tempfile.TemporaryDirectory(prefix='toe-public-replay-') as d:
        work=Path(d)
        shutil.copytree(ROOT/'baseline',work/'baseline')
        shutil.copytree(ROOT/'checkpoint_2026_09_06',work/'checkpoint_2026_09_06')
        shutil.copytree(ROOT/'checkpoint_2026_10_01',work/'checkpoint_2026_10_01')
        shutil.copytree(ROOT/'checkpoint_2026_10_15',work/'checkpoint_2026_10_15')
        shutil.copytree(ROOT/'checkpoint_2026_10_16',work/'checkpoint_2026_10_16')
        baseline=Path('baseline/TOE_N00AK_r1_PUBLIC')
        code=baseline/'code'
        cp=Path('checkpoint_2026_09_06')
        stages=[
            ('baseline input manifest',baseline/'verify_manifest.py',[]),
            ('baseline tree and conditional algebra',code/'verify_public.py',['--output','results/replay_normal.json']),
            ('baseline optimized semantic comparison',code/'verify_public.py',['--output','results/replay_optimized.json','--compare','results/replay_normal.json']),
            ('routed seed exact rational controls',cp/'routed_loop/verify_routed_pole.py',[]),
            ('separate tensor routing',cp/'routed_loop_independent.py',[]),
            ('exact weak-pulse lower bound',cp/'dynamic_eft_exact_bound.py',[]),
            ('scalar continuum controls',cp/'continuum_bridge/verify_continuum_bridge.py',[]),
            ('local Green representative and gauge ambiguity',Path('checkpoint_2026_10_01/verify_green_representative.py'),[]),
            ('sourced sterile contacts and Schur invariance',Path('checkpoint_2026_10_01/verify_sourced_sterile_identity.py'),[]),
            ('one-hypercharge-leg seed and transverse residue',Path('checkpoint_2026_10_15/verify_gauge_leg.py'),[]),
            ('independent gauge-leg angular reduction',Path('checkpoint_2026_10_15/verify_gauge_asymptotic.py'),[]),
            ('literal action and charge-conjugated gauge sign',Path('checkpoint_2026_10_15/operator_map/verify_action_sign.py'),[]),
            ('local transverse operators and sourced contacts',Path('checkpoint_2026_10_15/operator_map/verify_operator_map.py'),[]),
            ('exact transverse Ward-kernel classification',Path('checkpoint_2026_10_15/operator_map/verify_transverse_space.py'),[]),
            ('seven two-hypercharge-leg simplex residues',Path('checkpoint_2026_10_15/two_gauge_leg/verify_two_b.py'),[]),
            ('independent two-hypercharge-leg angular residues',Path('checkpoint_2026_10_15/two_gauge_leg/verify_two_b_angular.py'),[]),
            ('unintegrated two-hypercharge-leg vertices and Ward identity',Path('checkpoint_2026_10_15/two_gauge_leg/verify_two_b_vertices.py'),[]),
            ('one-weak-field index contractions and triangle residues',Path('checkpoint_2026_10_15/electroweak_source/verify_one_w.py'),[]),
            ('independent weak-field angular coefficients and photon projection',Path('checkpoint_2026_10_15/electroweak_source/verify_su2_independent.py'),[]),
            ('combined hypercharge source packet and equation residuals',Path('checkpoint_2026_10_15/electroweak_source/verify_selected_packet.py'),[]),
            ('indexed field redefinitions and integrated residual adjoints',Path('checkpoint_2026_10_15/field_redefinition/verify_field_redefinition.py'),[]),
            ('independent Grassmann and action-variation audit',Path('checkpoint_2026_10_15/field_redefinition/independent_action_audit.py'),[]),
            ('neutral quadratic projection and canonical pole/source transport',Path('checkpoint_2026_10_15/field_redefinition/verify_canonical.py'),[]),
            ('dimension-six normalization acting on a dimension-five source',Path('checkpoint_2026_10_15/field_redefinition/verify_lower_order_closure.py'),[]),
            ('independent component Grassmann lower-order closure',Path('checkpoint_2026_10_15/field_redefinition/independent_grassmann_closure_check.py'),[]),
            ('ordered two-weak and mixed-gauge simplex residues',Path('checkpoint_2026_10_15/nonabelian_gauge/verify_two_weak.py'),['--json']),
            ('independent non-Abelian two-gauge angular residues',Path('checkpoint_2026_10_15/nonabelian_gauge/verify_two_gauge_angular.py'),[]),
            ('matrix curvature and two/three-gauge functional Ward identities',Path('checkpoint_2026_10_15/nonabelian_gauge/verify_covariant_action.py'),[]),
            ('thirty-three three-gauge loop residues and neutral charged projection',Path('checkpoint_2026_10_15/nonabelian_gauge/verify_three_gauge.py'),[]),
            ('independent three-gauge exact cubature integration',Path('checkpoint_2026_10_15/nonabelian_gauge/verify_three_gauge_cubature.py'),[]),
            ('individual three-gauge diagram cross-route comparison',Path('checkpoint_2026_10_15/nonabelian_gauge/compare_three_gauge_routes.py'),[]),
            ('matrix electroweak source reduction and indexed field changes',Path('checkpoint_2026_10_15/interacting_reduction/verify_electroweak_reduction.py'),[]),
            ('independent non-Abelian action and combined weak-flavor audit',Path('checkpoint_2026_10_15/interacting_reduction/independent_ew_audit.py'),[]),
            ('ordered weak-current descendants and gauge-coordinate projections',Path('checkpoint_2026_10_15/interacting_reduction/verify_current_projection.py'),[]),
            ('interacting differential normalization and mass-source transport',Path('checkpoint_2026_10_15/interacting_reduction/verify_interacting_canonical.py'),[]),
            ('active and retained-mass strength algebraic bases',Path('checkpoint_2026_10_15/operator_basis/verify_active_tensor_basis.py'),[]),
            ('independent active tensor invariant and Grassmann quotient',Path('checkpoint_2026_10_15/operator_basis/independent_tensor_basis.py'),[]),
            ('spectator-current algebraic bases and identical-lepton Fierz transport',Path('checkpoint_2026_10_15/operator_basis/verify_current_operator_basis.py'),[]),
            ('independent current coefficient reconstruction and rank witnesses',Path('checkpoint_2026_10_15/operator_basis/independent_current_basis.py'),[]),
            ('symmetric-seed source image and relative flavor bound',Path('checkpoint_2026_10_15/operator_basis/verify_seed_flavor_image.py'),[]),
            ('selected all-field bridge and off-shell scalar Hessian',Path('checkpoint_2026_10_16/action_completion/verify_field_map_composition.py'),[]),
            ('independent neutral Grassmann action-map composition',Path('checkpoint_2026_10_16/action_completion/independent_composition_audit.py'),[]),
            ('scalar gauge-current reduction with interaction sources',Path('checkpoint_2026_10_16/action_completion/verify_scalar_current_reduction.py'),[]),
            ('independent scalar invariant and sourced-divergence audit',Path('checkpoint_2026_10_16/action_completion/independent_scalar_audit.py'),[]),
            ('charged-lepton and quark spectator Yukawa transport',Path('checkpoint_2026_10_16/action_completion/verify_spectator_yukawa.py'),[]),
            ('curvature scalar-derivative basis and row-source transport',Path('checkpoint_2026_10_16/derivative_projection/verify_curvature_derivative.py'),[]),
            ('independent curvature invariant and source audit',Path('checkpoint_2026_10_16/derivative_projection/independent_curvature_audit.py'),[]),
            ('complete neutrino-source active field-content projection',Path('checkpoint_2026_10_16/derivative_projection/verify_active_source_projection.py'),[]),
            ('independent active projection and sourced noncommutation audit',Path('checkpoint_2026_10_16/derivative_projection/independent_active_audit.py'),[]),
            ('complex Majorana mass recurrences and source covariance',Path('checkpoint_2026_10_16/derivative_projection/verify_complex_mass_covariance.py'),[]),
        ]
        if args.full:
            stages.extend([
                ('six pulse cases',cp/'fermion_portal/benchmark.py',[]),
                ('separate rotated-basis Radau review',cp/'fermion_portal/independent_review.py',[]),
            ])
        for name,script,extra in stages:
            flags=['-B']
            if name=='baseline optimized semantic comparison':flags+=['-O']
            command=[sys.executable,*flags,str(work/script),*extra]
            run=subprocess.run(command,cwd=(work/script).parent,text=True,capture_output=True,timeout=1200)
            result={'name':name,'script':script.as_posix(),'python_flags':flags,'exit_code':run.returncode,'stdout':run.stdout.strip()}
            results.append(result)
            print(name+': '+('PASS' if run.returncode==0 else 'FAIL'),flush=True)
            if run.returncode:
                print(run.stderr,file=sys.stderr)
                raise SystemExit(run.returncode)
        import numpy, sympy, mpmath
        environment={'python':platform.python_version(),'numpy':numpy.__version__,'sympy':sympy.__version__,'mpmath':mpmath.__version__,'platform':platform.system()+' '+platform.machine()}
        if args.full:
            import scipy
            environment['scipy']=scipy.__version__
        numerical={}
        if args.full:
            for key,rel in [('benchmark','fermion_portal/BENCHMARK_RESULTS.json'),('separate_review','fermion_portal/INDEPENDENT_REVIEW_RESULTS.json')]:
                numerical[key]=json.loads((work/cp/rel).read_text())
        release=json.loads((ROOT/'PUBLIC_MANIFEST.json').read_text())['release']
        report={'release':release,'status':'PASS','utc_completed':datetime.now(timezone.utc).isoformat(),'mode':'full' if args.full else 'quick','environment':environment,'payload_files_verified':count,'stages':results,'scientific_scope':'Declared tree/conditional algebra, one seed polynomial, local operator representative, sourced-contact and field-redefinition identities, one-hypercharge-leg transverse UV residue, literal action signs, restricted transverse tensor space and sourced flavor/current identities, seven-diagram two-hypercharge-leg UV pole matched to that action, unintegrated Ward controls, one-W weak-index and triangle matching, restricted neutral-photon pole cancellation, combined hypercharge source packet with off-shell equation-residual certificate, explicit selected first-order field substitutions and independent Grassmann action audit, original neutral quadratic projection with canonical pole/source transport, declared dimension-six times dimension-five contact benchmark and independent component audit, ordered two-W/mixed-BW pole integration with independent angular checks, thirty-three integrated three-gauge matter diagrams with independent exact cubature and per-diagram cross-route comparison, matrix-curvature functional Ward checks, restricted gauge-multiplicity degree closure and local neutral charged-weak projection, selected matrix electroweak source reduction with indexed real weak-field changes, combined weak-flavor Grassmann projection and current descendants, local differential normalization of the actual packet with explicit interaction remainder and differential mass-source transport, independent algebraic active and retained-mass strength sectors with ordered source transport, spectator-current Lorentz and weak invariant bases with identical-lepton Fierz flavor routing, separate current component-rank and formal source reconstruction, symmetric-seed flavor-image constraints and sharp relative algebraic source norm bound, selected all-field source-to-normalizing bridge with off-shell scalar-tadpole Hessian, independent neutral Grassmann full-action comparison, scalar gauge-current IBP and sourced-EOM transport, independent scalar invariant/component audit, charged-lepton and quark spectator Yukawa substitutions, curvature derivative-placement invariant spaces and source-preserving row transport, complete declared neutrino-source active field-content projection with independent Fourier/Grassmann audit, noncommuting projection and scalar-current EOM map, complex-Majorana sourced recurrences and flavor-coordinate covariance, inherited scalar comparisons, prescribed-background dynamics. No finite-photon prediction, complete physical matching, independent human peer review or empirical validation.','publication_replay_assertions_active':True,'temporary_copy_used':True,'numerical_results':numerical}
        if args.output:
            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(f'All {len(results)+1} stages passed; shipped files unchanged.',flush=True)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
