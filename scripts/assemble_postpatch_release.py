#!/usr/bin/env python3
"""Assemble the canonical post-patch scientific reproducibility release.

The release contains scientific code, data, numerical outputs, scientific audit
records, and reproducible figures/figure code. Manuscripts and all LaTeX build
artifacts are forbidden.
"""
from __future__ import annotations
import argparse, hashlib, json, shutil, tarfile, tempfile
from pathlib import Path

FORBIDDEN_EXT={'.tex','.bib','.bbl','.bcf','.blg','.aux','.toc','.lof','.lot','.run.xml','.synctex.gz','.fdb_latexmk','.fls'}
FORBIDDEN_PDF_PREFIX=('main','manuscript','paper','supplement')
FORBIDDEN_PATH_PARTS={'build_logs','clean_run_logs','latex','manuscript','paper_source','paper_sources'}
FORBIDDEN_BASENAMES={
    'verify_manuscript_data.py','audit_repository.py','generate_freeflux_runtime_table.py',
    'build_all.sh','reproduce_all.sh','run_full_audit.sh','run_full_audit_picasso.slurm'
}

def sha256(path: Path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()

def is_forbidden(path: Path):
    name=path.name.lower()
    parts={x.lower() for x in path.parts}
    if path.suffix.lower() in FORBIDDEN_EXT: return True
    if name in FORBIDDEN_BASENAMES: return True
    if parts & FORBIDDEN_PATH_PARTS: return True
    if path.suffix.lower()=='.pdf' and name.startswith(FORBIDDEN_PDF_PREFIX): return True
    if name.startswith('main_') and path.suffix.lower()=='.log': return True
    if name=='biber.log': return True
    return False

def purge_forbidden(root: Path):
    removed=[]
    for p in sorted(root.rglob('*'), reverse=True):
        if p.is_file() and is_forbidden(p):
            removed.append(str(p.relative_to(root))); p.unlink()
    for p in sorted(root.rglob('*'), reverse=True):
        if p.is_dir() and not any(p.iterdir()): p.rmdir()
    return sorted(removed)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--snapshot',required=True)
    ap.add_argument('--output-dir',required=True)
    ap.add_argument('--name',default='asymptotic_identifiability_reproducibility_v1.0.0')
    a=ap.parse_args()
    src=Path(__file__).resolve().parents[1]
    snap=Path(a.snapshot).resolve(); out=Path(a.output_dir).resolve(); out.mkdir(parents=True,exist_ok=True)
    if not snap.is_file(): raise SystemExit(f'snapshot not found: {snap}')
    with tempfile.TemporaryDirectory() as td0:
        td=Path(td0); rel=td/a.name
        shutil.copytree(src,rel,ignore=shutil.ignore_patterns('__pycache__','*.pyc','audit/full_runs'))
        ext=td/'snapshot'; ext.mkdir()
        with tarfile.open(snap,'r:gz') as tf:
            for m in tf.getmembers():
                p=Path(m.name)
                if p.is_absolute() or '..' in p.parts: raise SystemExit(f'unsafe tar member: {m.name}')
            tf.extractall(ext)
        # Only scientific trees are overlaid from the successful final snapshot.
        for d in ('validation','figures'):
            sp=ext/d
            if sp.exists(): shutil.copytree(sp,rel/d,dirs_exist_ok=True)
        runs=sorted((ext/'audit/full_runs').glob('*')) if (ext/'audit/full_runs').exists() else []
        if not runs: raise SystemExit('snapshot has no audit/full_runs entry')
        final_run=runs[-1]
        dst=rel/'audit'/'final_run_2279016'; shutil.copytree(final_run,dst,dirs_exist_ok=True)
        rr=dst/'science_run_result.json'
        if not rr.exists(): raise SystemExit('science_run_result.json missing')
        ev=json.loads(rr.read_text())
        if ev.get('exit_code')!=0 or not ev.get('completion_marker_found'):
            raise SystemExit(f'final run is not successful: {ev}')
        removed=purge_forbidden(rel)
        bad=[str(p.relative_to(rel)) for p in rel.rglob('*') if p.is_file() and is_forbidden(p)]
        if bad: raise SystemExit('forbidden paper/build artifacts remain: '+repr(bad))
        files=[]
        for p in sorted(x for x in rel.rglob('*') if x.is_file()):
            if p.name in {'MANIFEST_SHA256.txt','PROVENANCE.json'}: continue
            files.append((str(p.relative_to(rel)),sha256(p)))
        (rel/'MANIFEST_SHA256.txt').write_text(''.join(f'{h}  {n}\n' for n,h in files))
        prov={
          'canonical_scientific_run':'2279016',
          'run_result':ev,
          'snapshot_sha256':sha256(snap),
          'file_count_excluding_manifest':len(files),
          'forbidden_artifacts_removed':removed,
          'policy':'scientific code + data + numerical outputs + scientific audit + figures/figure code only; no manuscript, bibliography, LaTeX sources, or LaTeX build artifacts',
        }
        (rel/'PROVENANCE.json').write_text(json.dumps(prov,indent=2)+'\n')
        archive=out/f'{a.name}.tar.gz'
        with tarfile.open(archive,'w:gz') as tf: tf.add(rel,arcname=a.name)
        digest=sha256(archive)
        (out/f'{archive.name}.sha256').write_text(f'{digest}  {archive.name}\n')
        print(f'RELEASE={archive}')
        print(f'SHA256={digest}')
        print(f'FILES={len(files)}')
        print('FINAL_RUN=2279016 VERIFIED')
        print('FORBIDDEN_PAPER_OR_LATEX_ARTIFACTS=0')
if __name__=='__main__': main()
