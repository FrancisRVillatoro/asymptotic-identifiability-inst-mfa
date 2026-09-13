# Block 5 — Experimental interpretation and feasibility audit

## Purpose

This block separates mathematical order lifting from experimental deployability. The previous draft contained two statements that were too strong: (i) direct oxaloacetate (OAA) measurement was treated nearly as impossible, and (ii) slow fumarate was described as the "cheapest" intervention without an explicit cost model. Both are corrected.

## Evidence audit

### OAA
Zimmermann, Sauer & Zamboni (Analytical Chemistry 2014, 86:3232–3237, DOI 10.1021/ac500472c) developed quench-coupled derivatization for alpha-keto acids and demonstrated quantitative OAA analysis and stable-isotope mass-isotopomer profiling. Direct OAA is therefore **specialized but feasible**, not impossible. It remains a poor first-line design in the toy benchmark because it requires dedicated chemistry while changing only the prefactor, not the asymptotic order.

### Fumarate
Absolute microbial metabolomics has directly quantified fumarate (Bennett et al., Nature Chemical Biology 2009, DOI 10.1038/nchembio.186), and LC-MS isotope-tracer workflows include TCA intermediates (Mackay et al., Methods in Enzymology 2015, DOI 10.1016/bs.mie.2015.05.016). Targeted Synechocystis LC-MS panels have also included fumarate (Shi et al., Frontiers in Microbiology 2017, DOI 10.3389/fmicb.2017.00280). Slow fumarate MID measurement is therefore the **highest-priority near-term dynamic augmentation** among those tested, conditional on analytical coverage. The manuscript now says "lowest incremental experimental disruption" rather than "cheapest".

### Millisecond Asp-fast experiment
Nöh et al. (Journal of Biotechnology 2007, DOI 10.1016/j.jbiotec.2006.11.015) demonstrated rapid microbial INST sampling with 20 samples over a 16-s transient and immediate methanol quenching. Xu et al. (Methods in Molecular Biology 2024, DOI 10.1007/978-1-0716-3802-6_17) describe 0.1–0.5 s rapid quenching for plant INST-13C experiments. The toy design requires 1.43–22.86 ms. Its earliest point is about 70x faster than 0.1 s and 350x faster than 0.5 s. It is therefore retained as a **conceptual controllability/observability demonstration**, not a routine protocol.

### Absolute pool measurements
Bennett et al. (Nature Protocols 2008, DOI 10.1038/nprot.2008.107) provide stable-isotope ratio-based absolute intracellular metabolite quantification, while Bennett et al. 2009 demonstrate broad absolute microbial metabolomics. These support the practical plausibility of calibrated absolute pool constraints. They do **not** establish the 10% log-error used in the finite-noise section for this exact benchmark; that value is now explicitly labeled a working precision target.

### Synechocystis F6P versus GAP
Shi et al. 2017 report both F6P and GAP in a targeted Synechocystis central-carbon LC-MS panel, supporting analytical detectability in the organism. However, Nam et al. (Journal of Chromatography A 2021, 1656:462531, DOI 10.1016/j.chroma.2021.462531) found glyceraldehyde 3-phosphate and DHAP unstable during freeze-thaw and long-term storage because of reversible isomerization. Thus F6P is classified as the more straightforward absolute-pool target; GAP requires dedicated stability/recovery validation. The theory does not privilege GAP: another independent observable spanning the same weak direction would serve the same order-lifting role.

## Final classification

1. **Near-term:** fumarate MID at existing time points.
2. **Plausible algebraic additions:** absolute AKG/Cit and F6P, with platform-specific calibration.
3. **Specialized:** direct OAA isotope analysis and absolute GAP.
4. **Conceptual / technology-development:** 1.43–22.86 ms Asp-fast sampling.

No universal monetary ranking is claimed.
