#!/usr/bin/env python3
from pathlib import Path
import hashlib,sys
ROOT=Path(__file__).resolve().parents[1]
ref=ROOT/'validation/v101_reference/legacy_v100_expected_sha256.txt'
bad=[]
for line in ref.read_text().splitlines():
    if not line.strip(): continue
    expected,rel=line.split(None,1); rel=rel.strip(); p=ROOT/rel
    if not p.exists(): bad.append((rel,'MISSING',expected)); continue
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    if h!=expected: bad.append((rel,h,expected))
    else: print(rel,'OK')
if bad:
    print('INHERITED_V100_HASH_CHECK=FAIL')
    for x in bad: print(x)
    raise SystemExit(1)
print('INHERITED_V100_HASH_CHECK=PASS')
