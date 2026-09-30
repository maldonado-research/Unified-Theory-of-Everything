#!/usr/bin/env python3
"""Check the listed public payload bytes. This does not validate scientific claims."""
from pathlib import Path, PurePosixPath
import hashlib, json
root = Path(__file__).resolve().parent
manifest = json.loads((root / 'PUBLIC_MANIFEST.json').read_text(encoding='utf-8'))
if manifest.get('schema') != 'toe-public-payload-sha256-v1':
    raise SystemExit('Unsupported manifest schema')
seen = set()
for row in manifest['files']:
    rel = row['path']
    p = PurePosixPath(rel)
    if p.is_absolute() or '..' in p.parts or rel in seen:
        raise SystemExit('Invalid or duplicate manifest path')
    seen.add(rel)
    target = root / rel
    if target.is_symlink() or not target.is_file():
        raise SystemExit('Missing or nonregular file: ' + rel)
    data = target.read_bytes()
    if len(data) != row['bytes'] or hashlib.sha256(data).hexdigest() != row['sha256']:
        raise SystemExit('Hash or size mismatch: ' + rel)
print(json.dumps({'payload_files_verified': len(seen), 'sha256_status': 'PASS',
                  'scientific_validity_inferred': False}))
