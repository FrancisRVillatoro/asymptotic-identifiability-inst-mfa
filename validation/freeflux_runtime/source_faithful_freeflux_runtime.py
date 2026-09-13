#!/usr/bin/env python3
"""
Source-faithful extraction of the FreeFlux 0.3.8 INST simulation path.

The implementation below transcribes the simulation-relevant classes and
algorithms from commit ec05c47bbc2e4ac58bb85d39408ff4ef4016a15c:
  core/{metabolite,reaction,emu,model}.py
  core/mdv.py
  analysis/{simulate,inst_simulate}.py
  utils/utils.py
  io/inputs.py

Only import/compatibility changes are made for the supported runtime
environments (including Picasso Python 3.9 / pandas 1.4 and newer stacks):
  * np.float -> float
  * scipy.linalg.pinv2 -> scipy.linalg.pinv
  * DataFrame.map -> map/applymap compatibility fallback
  * OS CPU-affinity request omitted
  * fitting/optimization/result-display layers omitted

The EMU decomposition, A/B/M matrix construction, natural-abundance model,
tracer MDV construction, initialization, and piecewise-linear INST propagator
are otherwise the corresponding FreeFlux algorithms.
"""
from __future__ import annotations

import re
from collections import ChainMap, Counter, OrderedDict, deque
from collections.abc import Iterable
from functools import lru_cache, reduce
from itertools import chain, combinations_with_replacement, product
from numbers import Real
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import expm, pinv
from scipy.special import comb
from sympy import Integer, Matrix, Symbol, lambdify, symbols

COMMIT = "ec05c47bbc2e4ac58bb85d39408ff4ef4016a15c"
VERSION = "0.3.8"


class Metabolite:
    def __init__(self, id, atoms=None):
        self.id = id
        self.atoms = atoms
        if self.atoms and isinstance(self.atoms, str):
            self.atoms = [self.atoms]
        self.host_reactions = None

    def __hash__(self):
        if self.atoms:
            return hash(self.id) + sum(hash(atoms) for atoms in self.atoms)
        return hash(self.id)

    def __eq__(self, other):
        if not isinstance(other, Metabolite):
            return False
        if self.atoms:
            return self.id == other.id and set(self.atoms) == set(other.atoms)
        return (not other.atoms) and self.id == other.id

    @property
    def atoms_info(self):
        if self.atoms:
            return {atoms: 1 / len(self.atoms) for atoms in self.atoms}
        return None

    @property
    def n_carbons(self):
        if self.atoms:
            return len(re.search(r"[a-z]+", self.atoms[0]).group())
        return None

    def __repr__(self):
        atom_str = "(" + ",".join(self.atoms) + ")" if self.atoms else ""
        return f"{self.__class__.__name__} {self.id}{atom_str}"


class EMU:
    def __init__(self, id, metabolite, atom_nos):
        self.id = id
        if isinstance(metabolite, Metabolite):
            self.metabolite = metabolite
            self.metabolite_id = metabolite.id
        elif isinstance(metabolite, str):
            self.metabolite = Metabolite(metabolite)
            self.metabolite_id = metabolite
        else:
            raise TypeError("metabolite must be Metabolite or str")
        if isinstance(atom_nos, list):
            self.atom_nos = sorted(atom_nos)
        elif isinstance(atom_nos, str):
            self.atom_nos = sorted(map(int, atom_nos))
        else:
            self.atom_nos = sorted(atom_nos)
        self.size = len(self.atom_nos)

    def __hash__(self):
        return hash(self.metabolite_id) + hash(sum(self.atom_nos))

    def __eq__(self, other):
        if isinstance(other, Iterable):
            return type(other)([self]) == other
        return (
            isinstance(other, EMU)
            and self.metabolite_id == other.metabolite_id
            and self.atom_nos == other.atom_nos
        )

    def __lt__(self, other):
        if isinstance(other, Iterable):
            return type(other)([self]) < other
        if self.metabolite_id != other.metabolite_id:
            return self.metabolite_id < other.metabolite_id
        return self.atom_nos < other.atom_nos

    def __gt__(self, other):
        if isinstance(other, Iterable):
            return type(other)([self]) > other
        if self.metabolite_id != other.metabolite_id:
            return self.metabolite_id > other.metabolite_id
        return self.atom_nos > other.atom_nos

    @property
    @lru_cache()
    def equivalent_atom_nos(self):
        if len(self.metabolite.atoms_info) == 1:
            return None
        ref_atoms, equiv_atoms = self.metabolite.atoms_info
        mapping = dict(zip(ref_atoms, range(1, len(ref_atoms) + 1)))
        equiv = sorted(mapping[equiv_atoms[no - 1]] for no in self.atom_nos)
        return None if equiv == self.atom_nos else equiv

    @property
    @lru_cache()
    def equivalent(self):
        equiv = self.equivalent_atom_nos
        if equiv:
            id = self.metabolite_id + "_" + "".join(map(str, equiv))
            return EMU(id, self.metabolite, equiv)
        return None

    def __repr__(self):
        return f"{self.__class__.__name__} {self.metabolite_id}_{''.join(map(str, self.atom_nos))}"


class Reaction:
    def __init__(self, id, reversible=True):
        self.id = id
        self.reversible = reversible
        self.substrates_info = pd.DataFrame(columns=["metab", "stoy"])
        self.products_info = pd.DataFrame(columns=["metab", "stoy"])
        if reversible:
            self.fflux = Symbol(id + "_f")
            self.bflux = Symbol(id + "_b")
        else:
            self.flux = Symbol(id)
        self.host_models = None

    def add_substrates(self, substrates, stoichiometry):
        if not isinstance(substrates, list):
            substrates, stoichiometry = [substrates], [stoichiometry]
        new = pd.DataFrame(
            {"metab": substrates, "stoy": np.asarray(stoichiometry, dtype=float)},
            index=[sub.id for sub in substrates],
        )
        self.substrates_info = pd.concat((self.substrates_info, new))
        for sub in substrates:
            if sub.host_reactions is None:
                sub.host_reactions = {self}
            else:
                sub.host_reactions.add(self)

    def add_products(self, products, stoichiometry):
        if not isinstance(products, list):
            products, stoichiometry = [products], [stoichiometry]
        new = pd.DataFrame(
            {"metab": products, "stoy": np.asarray(stoichiometry, dtype=float)},
            index=[pro.id for pro in products],
        )
        self.products_info = pd.concat((self.products_info, new))
        for pro in products:
            if pro.host_reactions is None:
                pro.host_reactions = {self}
            else:
                pro.host_reactions.add(self)

    @property
    @lru_cache()
    def substrates(self):
        return sorted(self.substrates_info.index.unique().tolist())

    @property
    @lru_cache()
    def products(self):
        return sorted(self.products_info.index.unique().tolist())

    @property
    @lru_cache()
    def substrates_with_atoms(self):
        return sorted({sub for sub, row in self.substrates_info.iterrows() if row["metab"].atoms})

    @property
    @lru_cache()
    def products_with_atoms(self):
        return sorted({pro for pro, row in self.products_info.iterrows() if row["metab"].atoms})

    def _atom_mapping(self, reactant):
        if reactant == "substrate":
            info = self.substrates_info.loc[self.substrates_with_atoms, "metab"]
        elif reactant == "product":
            info = self.products_info.loc[self.products_with_atoms, "metab"]
        else:
            raise ValueError
        atom_info_all = []
        for _, metab in info.items():
            scenarios = []
            for atoms, coe in metab.atoms_info.items():
                scenarios.append({atom: [metab, no + 1, coe] for no, atom in enumerate(atoms)})
            atom_info_all.append(scenarios)
        raw = list(product(*atom_info_all))
        return [ChainMap(*scenario) for scenario in raw]

    @property
    @lru_cache()
    def _substrates_atom_mapping(self):
        return self._atom_mapping("substrate") if self.substrates_with_atoms else None

    @property
    @lru_cache()
    def _products_atom_mapping(self):
        return self._atom_mapping("product") if self.products_with_atoms else None

    def _find_precursor_EMUs(self, emu, direction="forward"):
        if self.reversible:
            if direction == "forward":
                atom_mapping = self._substrates_atom_mapping
            elif direction == "backward":
                atom_mapping = self._products_atom_mapping
            else:
                raise ValueError
        else:
            if direction != "forward":
                raise ValueError
            atom_mapping = self._substrates_atom_mapping

        raw = []
        for scenario in atom_mapping:
            for atoms, coe in emu.metabolite.atoms_info.items():
                pre_atoms = {}
                uni_coe = coe
                for atom in [atoms[no - 1] for no in emu.atom_nos]:
                    pre, pre_no, pre_coe = scenario[atom]
                    if pre not in pre_atoms:
                        uni_coe *= pre_coe
                        pre_atoms[pre] = [pre_no]
                    else:
                        pre_atoms[pre].append(pre_no)
                pre_emus = [
                    EMU(pre.id + "_" + "".join(map(str, sorted(nos))), pre, nos)
                    for pre, nos in pre_atoms.items()
                ]
                raw.append([pre_emus, uni_coe])
        counters = [Counter({tuple(sorted(pre)): coe}) for pre, coe in raw]
        merged = reduce(lambda x, y: x + y, counters)
        return [[list(pre), coe] for pre, coe in merged.items()]

    def __repr__(self):
        arrow = "<->" if self.reversible else "->"
        return f"Reaction {self.id}: {'+'.join(self.substrates)}{arrow}{'+'.join(self.products)}"


natAbuns = {
    "H": [0.999885, 0.000115],
    "C": [0.9893, 0.0107],
    "N": [0.99632, 0.00368],
    "O": [0.99757, 0.00038, 0.00205],
    "Si": [0.922297, 0.046832, 0.030872],
    "S": [0.9493, 0.0076, 0.0429, 0.0002],
}


class MDV:
    def __init__(self, fractions, nonnegative=True, normalize=True, base_atom="C"):
        self.value = np.asarray(fractions, dtype=float).copy()
        if nonnegative:
            self.value[self.value < 0] = 0
        if normalize and self.value.sum() != 0:
            self.value /= self.value.sum()
        self.base_atom = base_atom

    def __iter__(self):
        return iter(self.value)

    def __getitem__(self, key):
        return self.value[key]

    def __array__(self, dtype=None, copy=None):
        arr = self.value if dtype is None else self.value.astype(dtype)
        return arr.copy() if copy else arr

    def conv(self, mdv):
        if not isinstance(mdv, MDV):
            mdv = MDV(mdv)
        return MDV(gen_conv(self.value, mdv.value), base_atom=self.base_atom)

    def __mul__(self, other):
        if isinstance(other, Real):
            return MDV(other * self.value, nonnegative=False, normalize=False, base_atom=self.base_atom)
        if isinstance(other, Iterable):
            return self.conv(other)
        return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __add__(self, other):
        if isinstance(other, MDV):
            return MDV(self.value + other.value, nonnegative=False, normalize=False, base_atom=self.base_atom)
        return NotImplemented

    def __radd__(self, other):
        if other == 0:
            return self
        return self + other

    @property
    @lru_cache()
    def n_atoms(self):
        return self.value.size - 1

    def __repr__(self):
        return f"MDV([{', '.join(map(str, self.value.round(6)))}])"


def _isotopomer_combination(n_atoms, n_natural_isotops):
    all_combos = combinations_with_replacement(range(n_natural_isotops), n_atoms)
    combos1 = {}
    for combo in all_combos:
        combos1.setdefault(sum(combo), []).append(combo)
    combos2 = OrderedDict()
    for mass in sorted(combos1):
        for combo in combos1[mass]:
            combos2.setdefault(mass, []).append(Counter(combo))
    return combos2


def get_natural_MDV(n_atoms, base_atom="C"):
    nat = natAbuns[base_atom]
    combos = _isotopomer_combination(n_atoms, len(nat))
    mdv = []
    for counters in combos.values():
        ele = 0.0
        for counter in counters:
            item = 1.0
            left = n_atoms
            for isotop, count in counter.items():
                item *= comb(left, count) * nat[isotop] ** count
                left -= count
            ele += item
        mdv.append(ele)
    return MDV(mdv, base_atom=base_atom)


def get_substrate_MDV(atom_nos, labeling_pattern, percentage, purity, label_atom="C"):
    if not isinstance(labeling_pattern, list):
        if re.match(r"^0+$", labeling_pattern):
            raise ValueError("use natural MDV for natural substrate")
        labeling_pattern = [labeling_pattern]
    if not isinstance(percentage, list):
        percentage = [percentage]
    if sum(percentage) > 1:
        raise ValueError
    if not isinstance(purity, list):
        purity = [purity]
    n_atoms = len(atom_nos)
    single = {}
    for pat, pur in zip(labeling_pattern, purity):
        single[pat] = {
            "1": MDV([1 - pur, pur], base_atom=label_atom),
            "0": get_natural_MDV(1, base_atom=label_atom),
        }
    total = MDV(np.zeros(n_atoms + 1), base_atom=label_atom)
    for pat, per, pur in zip(labeling_pattern, percentage, purity):
        mdv = MDV([1], base_atom=label_atom)
        for atom_no in atom_nos:
            mdv *= single[pat][pat[atom_no - 1]]
        total += per * mdv
    return total + (1 - sum(percentage)) * get_natural_MDV(n_atoms, base_atom=label_atom)


def gen_conv(arr1, arr2):
    n1, n2 = len(arr1) - 1, len(arr2) - 1
    if n2 > n1:
        arr1, arr2, n1, n2 = arr2, arr1, n2, n1
    out = []
    for i in range(n1 + n2 + 1):
        if i <= n2:
            out.append(sum(arr1[i - j] * arr2[j] for j in range(i + 1)))
        elif i <= n1:
            out.append(sum(arr1[i - j] * arr2[j] for j in range(n2 + 1)))
        else:
            out.append(sum(arr1[i - j] * arr2[j] for j in range(i - n1, n2 + 1)))
    return np.asarray(out)


def conv(mdv1, mdv2):
    if not isinstance(mdv1, MDV):
        mdv1 = MDV(mdv1)
    if not isinstance(mdv2, MDV):
        mdv2 = MDV(mdv2)
    return mdv1 * mdv2


class Model:
    def __init__(self, name="unnamed"):
        self.name = name
        self.reactions_info = OrderedDict()
        self.target_EMUs = []
        self.timepoints = []
        self.substrate_MDVs = {}
        self.EAMs = {}
        self.matrix_As = {}
        self.matrix_Bs = {}
        self.matrix_Ms = {}
        self.initial_matrix_Xs = {}
        self.initial_matrix_Ys = {}
        self.initial_sim_MDVs = {}
        self.label_atom = None
        self.labeling_strategy = {}
        self.total_fluxes = pd.Series(dtype=float)
        self.concentrations = pd.Series(dtype=float)

    def add_reactions(self, reactions):
        if not isinstance(reactions, list):
            reactions = [reactions]
        self.reactions_info.update(OrderedDict((rxn.id, rxn) for rxn in reactions))
        for rxn in reactions:
            if rxn.host_models is None:
                rxn.host_models = {self}
            else:
                rxn.host_models.add(self)

    def read_from_file(self, file):
        data_raw = pd.read_csv(
            file, sep="\t", comment="#", header=None,
            names=["subs", "pros", "rev"], index_col=0,
        ).dropna()
        data = pd.DataFrame()
        for col, ser in data_raw.items():
            data[col] = ser.str.replace(r"\s+", "", regex=True) if ser.dtype == object else ser
        pattern_re = re.compile(r"([0-9\.]+|)([.\w]+)\(?([a-z0-9\.,]+|)\)?")
        for rxn, (subs_str, pros_str, rev) in data.iterrows():
            v = Reaction(rxn, reversible=bool(int(rev)))
            for stoy, metab, atoms in pattern_re.findall(subs_str):
                v.add_substrates(Metabolite(metab, atoms.split(",") if atoms else None), float(stoy) if stoy else 1.0)
            for stoy, metab, atoms in pattern_re.findall(pros_str):
                v.add_products(Metabolite(metab, atoms.split(",") if atoms else None), float(stoy) if stoy else 1.0)
            self.add_reactions(v)

    @property
    def metabolites_info(self):
        out = {}
        for rxn in self.reactions_info.values():
            for subid, sub in rxn.substrates_info["metab"].items():
                out.setdefault(subid, []).append(sub)
            for proid, pro in rxn.products_info["metab"].items():
                out.setdefault(proid, []).append(pro)
        return {k: list(set(v)) for k, v in out.items()}

    @property
    def metabolites(self):
        ids = []
        for rxn in self.reactions_info.values():
            ids += rxn.substrates + rxn.products
        return sorted(set(ids))

    @property
    def metabolites_with_atoms(self):
        ids = []
        for rxn in self.reactions_info.values():
            ids += rxn.substrates_with_atoms + rxn.products_with_atoms
        return sorted(set(ids))

    @property
    def reactions(self):
        return list(self.reactions_info)

    @lru_cache()
    def _full_net_stoichiometric_matrix(self, metabolites, reactions):
        net = pd.DataFrame(0.0, index=metabolites, columns=reactions)
        for rxnid, rxn in self.reactions_info.items():
            for sub in rxn.substrates:
                value = rxn.substrates_info.loc[sub, "stoy"]
                net.loc[sub, rxnid] = -value.sum() if isinstance(value, pd.Series) else -value
            for pro in rxn.products:
                value = rxn.products_info.loc[pro, "stoy"]
                net.loc[pro, rxnid] = value.sum() if isinstance(value, pd.Series) else value
        return net

    @lru_cache()
    def _full_total_stoichiometric_matrix(self, metabolites, reactions):
        net = self._full_net_stoichiometric_matrix(metabolites, reactions)
        cols, names = [], []
        for rxn, col in net.items():
            if self.reactions_info[rxn].reversible:
                cols += [col, -col + 0.0]
                names += [rxn + "_f", rxn + "_b"]
            else:
                cols.append(col); names.append(rxn)
        return pd.DataFrame(cols, index=names).T

    @property
    def end_substrates(self):
        total = self._full_total_stoichiometric_matrix(tuple(self.metabolites), tuple(self.reactions))
        nneg = (total < 0).sum(axis=1); npos = (total > 0).sum(axis=1)
        return sorted(total.index[(nneg > 0) & (npos == 0)].tolist())

    @property
    def end_products(self):
        total = self._full_total_stoichiometric_matrix(tuple(self.metabolites), tuple(self.reactions))
        nneg = (total < 0).sum(axis=1); npos = (total > 0).sum(axis=1)
        return sorted(total.index[(npos > 0) & (nneg == 0)].tolist())

    @property
    @lru_cache()
    def totalfluxids(self):
        ids = []
        for rxnid, rxn in self.reactions_info.items():
            ids.extend([rxnid + "_f", rxnid + "_b"] if rxn.reversible else [rxnid])
        return ids

    @property
    @lru_cache()
    def metabolite_adjacency_matrix(self):
        M = pd.DataFrame(index=self.metabolites_with_atoms, columns=self.metabolites_with_atoms)
        # Compatibility with the pinned Picasso environment (pandas 1.4.4):
        # DataFrame.map is unavailable there; DataFrame.applymap performs the
        # same elementwise initialization. Prefer map on newer pandas.
        if hasattr(M, "map"):
            M = M.map(lambda _: [])
        else:
            M = M.applymap(lambda _: [])
        for rxn in self.reactions_info.values():
            for sub in rxn.substrates_with_atoms:
                for pro in rxn.products_with_atoms:
                    if rxn.reversible:
                        M.loc[sub, pro].append(rxn)
                        M.loc[pro, sub].append(rxn)
                    else:
                        M.loc[sub, pro].append(rxn)
        return M

    def _BFS(self, iniEMU):
        MAM = self.metabolite_adjacency_matrix
        info = {}
        searched = []
        to_search = deque([iniEMU])
        while to_search:
            current = to_search.pop()
            searched.append(current)
            forming = list(set(chain(*[cell for cell in MAM[current.metabolite_id] if cell])))
            for rxn in forming:
                if rxn.reversible:
                    if current.metabolite_id in rxn.products_with_atoms:
                        as_pro = rxn.products_info["metab"][current.metabolite_id]
                        direction, flux = "forward", rxn.fflux
                    else:
                        as_pro = rxn.substrates_info["metab"][current.metabolite_id]
                        direction, flux = "backward", rxn.bflux
                else:
                    as_pro = rxn.products_info["metab"][current.metabolite_id]
                    direction, flux = "forward", rxn.flux
                if isinstance(as_pro, pd.Series):
                    offset = 1 / as_pro.size
                    as_pro = list(as_pro)
                else:
                    offset = 1.0
                    as_pro = [as_pro]
                for metab in as_pro:
                    current_variant = EMU(current.id, metab, current.atom_nos)
                    for pre_emus, coe in rxn._find_precursor_EMUs(current_variant, direction):
                        for pre in pre_emus:
                            if pre not in searched and pre not in to_search:
                                to_search.appendleft(pre)
                        info.setdefault(current_variant.size, []).append(
                            [current_variant, pre_emus, offset * coe * flux]
                        )
        return info

    def _get_original_EAMs(self, iniEMU):
        info = self._BFS(iniEMU)
        out = {}
        for size, rows in info.items():
            non_sources = set(row[0] for row in rows)
            sources = sorted(
                set(tuple(row[1]) if len(row[1]) > 1 else row[1][0] for row in rows) - non_sources
            )
            EAM = pd.DataFrame(Integer(0), index=sorted(non_sources) + sources, columns=sorted(non_sources))
            for emu, pre, flux in rows:
                idx = pre[0] if len(pre) == 1 else tuple(pre)
                EAM.loc[[idx], emu] += flux
            out[size] = EAM
        return out

    def _combine_equivalent_EMUs(self, EAMs):
        out = {}
        for size, EAM in EAMs.items():
            combined = EAM.copy(deep=True)
            done = []
            for emu in list(combined.columns):
                if emu in done or emu not in combined.columns:
                    continue
                equiv = emu.equivalent
                if equiv in combined.columns:
                    combined.loc[:, emu] = combined.loc[:, [emu, equiv]].sum(axis=1) / 2
                    combined.drop(equiv, axis=1, inplace=True)
                    combined.loc[emu, :] = combined.loc[[emu, equiv], :].sum()
                    combined.drop(equiv, inplace=True)
                    done.append(equiv)
            out[size] = combined
        return out

    def get_emu_adjacency_matrices(self, iniEMU, lump=True):
        # INST uses lump=False; only this branch is needed here.
        ori = self._get_original_EAMs(iniEMU)
        if lump:
            raise NotImplementedError("This runtime extraction is for INST (lump=False)")
        return self._combine_equivalent_EMUs(ori)

    def _merge_EAMs(self, EAM1, EAM2):
        non_sources = EAM2.columns.union(EAM1.columns)
        sources = EAM2.index.difference(EAM2.columns).union(EAM1.index.difference(EAM1.columns))
        merged = pd.DataFrame(Integer(0), index=non_sources.append(sources), columns=non_sources)
        merged.loc[EAM1.index, EAM1.columns] = EAM1
        merged.loc[EAM2.index, EAM2.columns] = EAM2
        return merged

    def _merge_all_EAMs(self, *all_eams):
        out = {}
        maxsize = max(max(eams) for eams in all_eams)
        for size in range(1, maxsize + 1):
            current = [eams.get(size) for eams in all_eams if isinstance(eams.get(size), pd.DataFrame)]
            if current:
                out[size] = reduce(self._merge_EAMs, current)
        return out

    def _decompose_network(self, metabolites, atom_nos, lump=True, n_jobs=1):
        emus = [EMU(m + "_" + nos, Metabolite(m), nos) for m, nos in zip(metabolites, atom_nos)]
        all_eams = [self.get_emu_adjacency_matrices(emu, lump) for emu in emus]
        return self._merge_all_EAMs(*all_eams)


class Calculator:
    def __init__(self, model):
        self.model = model

    def _calculate_substrate_MDVs(self, extra_subs=None):
        extra_subs = [] if extra_subs is None else list(extra_subs)
        for size in self.model.matrix_Bs:
            for source in self.model.matrix_Bs[size][2]:
                iterable = source if isinstance(source, Iterable) else (source,)
                for emu in iterable:
                    metabid = emu.metabolite_id
                    if metabid in self.model.end_substrates + extra_subs:
                        if metabid in self.model.labeling_strategy:
                            pat, per, pur = self.model.labeling_strategy[metabid]
                            self.model.substrate_MDVs[emu] = get_substrate_MDV(
                                emu.atom_nos, pat, per, pur, label_atom=self.model.label_atom
                            )
                        else:
                            self.model.substrate_MDVs[emu] = get_natural_MDV(emu.size, base_atom=self.model.label_atom)

    def _lambdify_matrix_As_and_Bs(self):
        for size, EAM in self.model.EAMs.items():
            pre = EAM.copy(deep=True)
            for emu in EAM.columns:
                pre.loc[emu, emu] = -pre[emu].sum()
            A = pre.loc[pre.columns, :].T
            B = -pre.loc[pre.index.difference(pre.columns), :].T
            matA, matB = Matrix(A), Matrix(B)
            fluxidsA = list(map(str, matA.free_symbols))
            fluxidsB = list(map(str, matB.free_symbols))
            self.model.matrix_As[size] = [lambdify(fluxidsA, matA, modules="numpy"), fluxidsA, A.columns.tolist()]
            self.model.matrix_Bs[size] = [lambdify(fluxidsB, matB, modules="numpy"), fluxidsB, B.columns.tolist()]

    def _lambdify_matrix_Ms(self):
        for size, EAM in self.model.EAMs.items():
            matM = Matrix(np.diag(symbols([emu.metabolite_id for emu in EAM.columns])))
            metabids = list(map(str, matM.free_symbols))
            self.model.matrix_Ms[size] = [lambdify(metabids, matM, modules="numpy"), metabids]

    def _calculate_initial_matrix_Xs(self):
        for size in self.model.matrix_As:
            n = len(self.model.matrix_As[size][2])
            self.model.initial_matrix_Xs[size] = np.vstack(
                [get_natural_MDV(size, base_atom=self.model.label_atom).value] * n
            )

    def _calculate_initial_matrix_Ys(self):
        for size in self.model.matrix_Bs:
            iniY = []
            for source in self.model.matrix_Bs[size][2]:
                if not isinstance(source, Iterable):
                    sourceMDV = self.model.substrate_MDVs[source]
                else:
                    mdvs = []
                    for emu in source:
                        mdvs.append(
                            get_natural_MDV(emu.size, base_atom=self.model.label_atom)
                            if emu not in self.model.substrate_MDVs
                            else self.model.substrate_MDVs[emu]
                        )
                    sourceMDV = reduce(conv, mdvs)
                iniY.append(np.asarray(sourceMDV))
            self.model.initial_matrix_Ys[size] = np.asarray(iniY)

    def _build_initial_sim_MDVs(self):
        for size in sorted(self.model.matrix_As):
            for emu, ini in zip(self.model.matrix_As[size][2], self.model.initial_matrix_Xs[size]):
                if emu.id in self.model.target_EMUs:
                    self.model.initial_sim_MDVs[emu.id] = {0: MDV(ini)}

    def _calculate_inst_MDVs(self):
        sim = {}
        Ys, Xs = {}, {}
        t1 = 0.0
        for size in sorted(self.model.matrix_As):
            Ys.setdefault(t1, {})[size] = self.model.initial_matrix_Ys[size]
            Xs.setdefault(t1, {})[size] = self.model.initial_matrix_Xs[size]
        for t in self.model.timepoints:
            if t == 0.0:
                continue
            t0, t1 = t1, t
            dt = t1 - t0
            for size in sorted(self.model.matrix_As):
                lambA, fluxidsA, productEMUs = self.model.matrix_As[size]
                lambB, fluxidsB, sourceEMUs = self.model.matrix_Bs[size]
                lambM, metabids = self.model.matrix_Ms[size]
                A = np.asarray(lambA(*self.model.total_fluxes[fluxidsA]), dtype=float)
                B = np.asarray(lambB(*self.model.total_fluxes[fluxidsB]), dtype=float)
                M = np.asarray(lambM(*self.model.concentrations[metabids]), dtype=float)
                Minv = pinv(M, check_finite=True)
                F = Minv @ A
                Finv = pinv(F, check_finite=True)
                I = np.eye(*F.shape)
                Phi = expm(F * dt)
                Gamma = (Phi - I) @ Finv
                Omega = (Gamma / dt - I) @ Finv
                X_t0 = Xs[t0][size]
                Y_t0 = Ys[t0][size]
                G_t0 = Minv @ B @ Y_t0
                Y_t1 = []
                for source in sourceEMUs:
                    if not isinstance(source, Iterable):
                        sourceMDV = self.model.substrate_MDVs[source]
                    else:
                        mdvs = []
                        for emu in source:
                            mdv = ChainMap(sim, self.model.substrate_MDVs)[emu]
                            if isinstance(mdv, dict):
                                mdv = mdv[t1]
                            mdvs.append(mdv)
                        sourceMDV = reduce(conv, mdvs)
                    Y_t1.append(np.asarray(sourceMDV))
                Y_t1 = np.asarray(Y_t1)
                G_t1 = Minv @ B @ Y_t1
                X_t1 = Phi @ X_t0 - Gamma @ G_t0 - Omega @ (G_t1 - G_t0)
                Ys.setdefault(t1, {})[size] = Y_t1
                Xs.setdefault(t1, {})[size] = X_t1
                for emu, mdv in zip(productEMUs, X_t1):
                    sim.setdefault(emu, {})[t1] = mdv
        return {emu.id: mdvs for emu, mdvs in sim.items()}


class InstSimulator:
    def __init__(self, model):
        self.model = model
        self.calculator = Calculator(model)

    def set_target_EMUs(self, target_emus):
        for metabid, atom_nos in target_emus.items():
            if isinstance(atom_nos, list):
                if any(isinstance(item, Iterable) for item in atom_nos):
                    for nos in atom_nos:
                        nos = nos if isinstance(nos, str) else "".join(map(str, nos))
                        self.model.target_EMUs.append(metabid + "_" + nos)
                else:
                    self.model.target_EMUs.append(metabid + "_" + "".join(map(str, atom_nos)))
            else:
                self.model.target_EMUs.append(metabid + "_" + atom_nos)

    def set_labeling_strategy(self, labeled_substrate, labeling_pattern, percentage, purity, label_atom="C"):
        self.model.labeling_strategy[labeled_substrate] = [labeling_pattern, percentage, purity]
        self.model.label_atom = label_atom

    def set_flux(self, fluxid, value):
        self.model.total_fluxes[fluxid] = value

    def set_concentration(self, metabid, value):
        self.model.concentrations[metabid] = value

    def set_timepoints(self, timepoints):
        self.model.timepoints = sorted(set(self.model.timepoints + list(timepoints)))
        if 0 not in self.model.timepoints:
            self.model.timepoints = [0] + self.model.timepoints

    def prepare(self):
        metabids, atom_nos = [], []
        for emuid in self.model.target_EMUs:
            metabid, nos = emuid.split("_")
            metabids.append(metabid); atom_nos.append(nos)
        self.model.EAMs = self.model._decompose_network(metabids, atom_nos, lump=False, n_jobs=1)
        self.calculator._lambdify_matrix_As_and_Bs()
        self.calculator._lambdify_matrix_Ms()
        self.calculator._calculate_substrate_MDVs()
        self.calculator._calculate_initial_matrix_Xs()
        self.calculator._calculate_initial_matrix_Ys()
        self.calculator._build_initial_sim_MDVs()

    def simulate_raw(self):
        raw = self.calculator._calculate_inst_MDVs()
        out = {}
        for emuid in self.model.target_EMUs:
            mdvs = {0.0: self.model.initial_sim_MDVs[emuid][0].value.copy()}
            mdvs.update({float(t): np.asarray(v, dtype=float) for t, v in raw[emuid].items()})
            out[emuid] = mdvs
        return out


def build_toy_model(reactions_file: str | Path):
    model = Model("toy")
    model.read_from_file(reactions_file)
    sim = InstSimulator(model)
    sim.set_target_EMUs({"Glu": [[1, 2, 3], "12345"], "Cit": "2345"})
    sim.set_labeling_strategy(
        "AcCoA", labeling_pattern=["01", "11"], percentage=[0.25, 0.25], purity=[1, 1], label_atom="C"
    )
    fluxes = {"v1": 10, "v2": 10, "v3": 5, "v4": 5, "v5": 5, "v6_f": 12.5, "v6_b": 7.5, "v7": 5}
    concentrations = {"OAA": 0.1, "Cit": 5, "AKG": 0.3, "Suc": 1, "Fum": 0.2, "Glu": 0.5}
    for k, v in fluxes.items(): sim.set_flux(k, v)
    for k, v in concentrations.items(): sim.set_concentration(k, v)
    return model, sim
