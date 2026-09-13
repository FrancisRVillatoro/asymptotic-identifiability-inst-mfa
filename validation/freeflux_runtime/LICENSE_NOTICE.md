# License and provenance notice

`source_faithful_freeflux_runtime.py` is a compatibility extraction derived
from the simulation-relevant source code of **FreeFlux 0.3.8**, authored by
Chao Wu and distributed under the GNU General Public License version 3.

Upstream source:

- repository: `https://github.com/Chaowu88/freeflux`
- commit: `ec05c47bbc2e4ac58bb85d39408ff4ef4016a15c`
- upstream license: GNU GPL v3.0

Compatibility changes made in the extraction:

1. deprecated `np.float` replaced by built-in `float`;
2. removed `scipy.linalg.pinv2` replaced by `scipy.linalg.pinv`;
3. fitting, optimization, plotting and result-display layers omitted because
   they do not enter the INST simulation calculation;
4. Linux CPU-affinity request omitted.

The EMU decomposition, equivalent-EMU combination, natural-abundance and
tracer-MDV calculations, construction of the A/B/M matrices, initialization,
and the piecewise-linear INST propagation algorithm are retained in
source-faithful form.

The complete GPL v3 license text is included as `COPYING.GPL-3`.
