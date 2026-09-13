# Block 3 — Automatic fast-block selection and partition robustness

## Goal

Replace manual choice of the fast block by a reproducible candidate-generation rule, then audit whether the asymptotic identifiability orders are stable to plausible changes of the cutoff.

## Automatic candidate rule

For positive turnover rates `k_i`, sort them increasingly and compute

`g_i = log10(k_(i+1)/k_i)`.

Only cuts that leave at most

`q_max = max(3, ceil(n/4))`

pools in the high-turnover suffix are considered. Among those, retain local maxima of `g_i`; the largest local maximum is the primary candidate. This is only a proposal stage. Candidate blocks are subsequently tested by (i) an EMU fast-block spectral diagnostic, (ii) an epsilon-homotopy singular-value sweep, (iii) output convergence to the reduced limit, and (iv) principal angles of weak subspaces under neighboring cutoffs.

Model-declared auxiliary/dummy atom-transfer pools are retained in the raw candidate but may be removed for a biochemical interpretation. In the Synechocystis benchmark, TA and TK are the auxiliary atom-transfer pools.

## Synechocystis turnover gaps

The top-tail local maxima are:

| cut_after   |   gap_log10 |   ratio |   raw_fast_size |   eligible_fast_size |
|:------------|------------:|--------:|----------------:|---------------------:|
| G3P         |   0.565587  | 3.67779 |               3 |                    2 |
| Mal         |   0.0924411 | 1.2372  |               7 |                    6 |
| TA          |   0.0226833 | 1.05362 |               1 |                    1 |

The primary gap is therefore the G3P -> F6P jump,

`ratio = 3.677789`,

which produces the raw fast block `{F6P, TA, GAP}`. After removal of the auxiliary TA pool, the biochemical fast block is `{F6P, GAP}`.

## Robustness sweep

The asymptotic sweep used epsilon = 1, 1/2, ..., 1/512 and fitted the last four points. Singular-value orders were classified as m=0 for |slope|<0.15, m=1 for 0.75<slope<1.25, and m=2 for 1.75<slope<2.25. One exact global flux--pool scale null was removed from the count.

| partition           | fast_pools                 |   n_fast |   m1_count | m1_exponents                                                   |   m2_count |   other_count |   pool_gap_ratio |   emu_to_slow_turnover_ratio |   output_convergence_order |   angle_to_primary_phys_max_deg | classification      |
|:--------------------|:---------------------------|---------:|-----------:|:---------------------------------------------------------------|-----------:|--------------:|-----------------:|-----------------------------:|---------------------------:|--------------------------------:|:--------------------|
| auto_secondary_phys | G2P,S7P,E4P,G3P,F6P,GAP    |        6 |          5 | 0.990065;0.999024;0.999846;0.998270;1.001231                   |          0 |             4 |         0.126946 |                     0.03027  |                   0.995808 |                     5.6254      | broad-cut/crossover |
| auto_secondary_raw  | G2P,S7P,E4P,G3P,F6P,TA,GAP |        7 |          7 | 0.782093;0.896114;0.993380;0.996735;1.002302;0.999715;1.004876 |          0 |             2 |         1.2372   |                     0.168711 |                   0.995609 |                     3.55296     | broad-cut/crossover |
| auto_fast1_raw      | GAP                        |        1 |          1 | 0.999517                                                       |          0 |             0 |         1.05362  |                     1.05362  |                   0.999575 |                     0.236015    | quantized-stable    |
| auto_primary_phys   | F6P,GAP                    |        2 |          2 | 0.999381;1.002753                                              |          0 |             0 |         0.977063 |                     0.13882  |                   0.998855 |                     1.91312e-14 | quantized-stable    |
| primary_plus1_phys  | G3P,F6P,GAP                |        3 |          3 | 0.999598;0.998760;1.002300                                     |          0 |             0 |         0.265666 |                     0.128688 |                   0.998123 |                     0.100739    | quantized-stable    |
| auto_primary_raw    | F6P,TA,GAP                 |        3 |          3 | 0.988309;0.999218;1.003367                                     |          0 |             0 |         3.67779  |                     0.42449  |                   0.998519 |                     1.29196     | quantized-stable    |
| primary_plus1_raw   | G3P,F6P,TA,GAP             |        4 |          4 | 0.987672;0.999609;0.997994;1.003510                            |          0 |             0 |         1.40563  |                     0.596675 |                   0.997806 |                     1.23741     | quantized-stable    |

### Primary-cut stability

The primary raw block gives three m=1 modes with exponents

`0.988309;0.999218;1.003367`.

Removing the auxiliary TA pool gives exactly two m=1 modes,

`0.999381;1.002753`,

and no m=2 mode. Adding the nearest slower biochemical pool G3P produces exactly one additional m=1 mode,

`0.999598;0.998760;1.002300`,

without materially rotating the original two-dimensional F6P/GAP weak subspace: the largest principal angle is only

`0.100739 degrees`.

For the raw block including TA, the corresponding angle is `1.291961 degrees`; for the one-boundary-expanded raw block it is `1.237407 degrees`.

Thus the physical F6P/GAP weak subspace is robust to the primary-cut perturbations; adding a pool adds a first-order mode rather than changing the order of the existing ones.

### Secondary broad cut is rejected as a clean two-scale description

The second local gap, after Mal, has only ratio `1.237203` and proposes seven raw fast pools. At epsilon down to 1/512 its spectrum is not cleanly quantized: the raw block has two crossover slopes outside the m=0,1,2 bands, and the biochemical six-pool block has four. Its principal-angle deviation from the primary biochemical weak subspace increases to 3.553 and 5.625 degrees, respectively.

The nominal EMU spectral audit also flags this broad cut: the ratio of the slowest restricted fast-EMU decay rate to the fastest remaining turnover is only 0.169 (raw) and 0.030 (biochemical). The primary raw candidate is also only an incipient, not fully separated, nominal fast block (`0.424`), which explains why the homotopy is needed to expose the asymptotic hierarchy.

This is an important distinction: the turnover-gap rule proposes a homotopy; it does not assert that epsilon=1 is already in a strongly separated asymptotic regime.

## Output reduction check

For all primary-cut variants, successive output differences converge with first-order slope very close to one:

- raw primary: 0.998519;
- biochemical primary: 0.998855;
- biochemical primary + G3P: 0.998123.

At the final halving step the relative output differences are O(1e-4) or smaller. This independently confirms that the same epsilon-homotopy which generates the information loss also approaches a well-defined reduced output.

## Cross-check on the FreeFlux toy TCA benchmark

Applying the same gap rule to the six published turnover rates gives:

| cut_after   |   gap_log10 |   ratio |   fast_size | fast_pools   |
|:------------|------------:|--------:|------------:|:-------------|
| Glu         |    0.522879 | 3.33333 |           3 | AKG,Fum,OAA  |
| Fum         |    0.447158 | 2.8     |           1 | OAA          |

The primary automatic cut is after Glu and returns

`{AKG, Fum, OAA}`,

exactly the three-pool fast block previously selected manually. For that benchmark the already-validated singular orders are `[0,0,0,1,1,2]` and the FIM exponents `[0,0,0,2,2,4]`. Hence the automatic candidate rule recovers the homotopy that exposes the two m=1 modes and the moment-degenerate m=2 OAA--Fum direction.

## Algorithmic conclusion

The practical rule emerging from the two benchmarks is:

1. compute physical turnover rates;
2. generate high-tail local-gap candidates automatically;
3. preserve a raw candidate and a semantically filtered biochemical candidate;
4. audit the actual EMU fast-block spectrum at the nominal point;
5. run an epsilon sweep and demand near-integer order stabilization;
6. perturb the cutoff by one boundary and compare principal angles of the weak subspace;
7. if neighboring cuts produce crossover slopes or large rotations, report the model as continuously/multiscale rather than force a binary fast--slow split.

The cutoff therefore becomes reproducible and falsifiable. The final asymptotic classification is not inferred from turnover gaps alone; it is accepted only after the singular-value and subspace-stability audits.

## Scientific interpretation

The Synechocystis primary block is robust in the sense relevant to the paper: F6P and GAP each generate a first-order information-loss direction, and nearby cutoff changes only append additional first-order modes. No m=2 mode appears. The broad secondary cut fails the clean two-scale audit. Together with the toy TCA result, this supports the interpretation that m=1 is the generic fast-turnover information loss, whereas m=2 requires an additional moment degeneracy rather than merely a larger fast block.
