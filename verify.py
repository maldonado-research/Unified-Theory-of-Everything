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
        import numpy
        environment={'python':platform.python_version(),'numpy':numpy.__version__,'platform':platform.system()+' '+platform.machine()}
        if args.full:
            import scipy
            environment['scipy']=scipy.__version__
        numerical={}
        if args.full:
            for key,rel in [('benchmark','fermion_portal/BENCHMARK_RESULTS.json'),('separate_review','fermion_portal/INDEPENDENT_REVIEW_RESULTS.json')]:
                numerical[key]=json.loads((work/cp/rel).read_text())
        report={'release':'TOE-GitHub-2026.09.30','status':'PASS','utc_completed':datetime.now(timezone.utc).isoformat(),'mode':'full' if args.full else 'quick','environment':environment,'payload_files_verified':count,'stages':results,'scientific_scope':'Declared tree/conditional algebra, one seed polynomial, inherited scalar comparisons, prescribed-background dynamics. No complete physical matching, independent human peer review or empirical validation.','publication_replay_assertions_active':True,'temporary_copy_used':True,'numerical_results':numerical}
        if args.output:
            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(f'All {len(results)+1} stages passed; shipped files unchanged.',flush=True)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
