#!/usr/bin/env python3
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
EXT={'.tex','.bib','.bbl','.bcf','.blg','.aux','.toc','.lof','.lot','.run.xml','.synctex.gz','.fdb_latexmk','.fls'}
BAD_NAMES={'verify_manuscript_data.py','audit_repository.py','generate_freeflux_runtime_table.py','build_all.sh','reproduce_all.sh','run_full_audit.sh','run_full_audit_picasso.slurm'}
BAD_PARTS={'build_logs','clean_run_logs'}
PDF_PREFIX=('main','manuscript','paper','supplement')
def bad(p):
    n=p.name.lower(); parts={x.lower() for x in p.parts}
    return (p.suffix.lower() in EXT or n in BAD_NAMES or bool(parts & BAD_PARTS) or
            (p.suffix.lower()=='.pdf' and n.startswith(PDF_PREFIX)) or
            (n.startswith('main_') and p.suffix.lower()=='.log') or n=='biber.log')
hits=[str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and bad(p)]
if hits:
    print('FORBIDDEN ARTIFACTS FOUND:')
    print('\n'.join(hits)); sys.exit(1)
print('NO_PAPER_ARTIFACTS=PASS')
