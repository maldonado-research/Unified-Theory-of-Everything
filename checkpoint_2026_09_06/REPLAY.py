"""Verify a clean copy of the focused TOE checkpoint and replay its checks."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

if sys.flags.optimize:
    raise SystemExit("Run REPLAY.py without -O so scientific assertions remain active.")

ROOT=Path(__file__).resolve().parent
manifest=json.loads((ROOT/'MANIFEST.json').read_text())
for item in manifest['files']:
    target=ROOT/item['path']
    if not target.is_file() or hashlib.sha256(target.read_bytes()).hexdigest()!=item['sha256']:
        raise SystemExit('Input checksum mismatch: '+item['path'])

stages=[
    ('routed UV polynomial','routed_loop/verify_routed_pole.py'),
    ('independent tensor routing','routed_loop_independent.py'),
    ('exact nonzero weak pulse bound','dynamic_eft_exact_bound.py'),
    ('continuum spectral checks','continuum_bridge/verify_continuum_bridge.py'),
    ('six fermion pulse cases','fermion_portal/benchmark.py'),
    ('independent rotated-basis solver','fermion_portal/independent_review.py'),
]
results=[]
for name,script in stages:
    run=subprocess.run([sys.executable,'-B',str(ROOT/script)],cwd=ROOT,
                       text=True,capture_output=True)
    results.append({'name':name,'script':script,'exit_code':run.returncode,
                    'stdout':run.stdout,'stderr':run.stderr})
    print(name+': '+('PASS' if run.returncode==0 else 'FAIL'),flush=True)
    if run.returncode:
        (ROOT/'REPLAY_RESULT.json').write_text(json.dumps({'status':'FAIL','stages':results},indent=2)+'\n')
        raise SystemExit(run.returncode)

# Recomputed numerical receipts may differ on other platforms. Record them;
# only scientific checks, not cross-platform bit identity, are required.
changed=[]
for item in manifest['files']:
    digest=hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()
    if digest!=item['sha256']:
        changed.append({'path':item['path'],'new_sha256':digest})
report={'status':'PASS','input_hashes_verified':len(manifest['files']),
        'verification_stages':len(results),'stages':results,
        'regenerated_files_with_changed_hashes':changed,
        'scope':'Bounded equations and toy models; not full dim7 matching, parent spectral proof, or empirical validation.'}
(ROOT/'REPLAY_RESULT.json').write_text(json.dumps(report,indent=2)+'\n')
print('All six stages passed.',flush=True)
