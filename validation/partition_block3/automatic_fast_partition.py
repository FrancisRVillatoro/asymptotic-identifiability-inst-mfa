#!/usr/bin/env python3
"""Automatic turnover-gap candidate generation for fast--slow INST-MFA homotopies.

The routine does not declare a partition 'correct'. It proposes high-turnover suffixes,
which must then pass a dynamical/spectral and asymptotic-order robustness audit.
"""
from __future__ import annotations
import math
from typing import Iterable
import numpy as np
import pandas as pd


def propose_fast_partitions(turnovers: pd.Series | dict,
                            max_fast_fraction: float = 0.25,
                            min_tail_size: int = 3,
                            auxiliary_pools: Iterable[str] = ()) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return sorted turnover table and top-tail local-gap candidate partitions.

    Parameters
    ----------
    turnovers
        Mapping/Series metabolite -> positive turnover rate.
    max_fast_fraction
        Only cuts leaving at most this fraction of the pools in the fast suffix are considered.
    min_tail_size
        Ensures small networks can still propose a three-pool fast block.
    auxiliary_pools
        Model-declared dummy/atom-transfer pools excluded only from the biochemical interpretation;
        raw candidates are always retained.
    """
    s=pd.Series(turnovers,dtype=float).sort_values()
    if (s<=0).any(): raise ValueError('turnovers must be positive')
    tab=pd.DataFrame({'metabolite':s.index.astype(str),'turnover':s.values})
    tab['log10_turnover']=np.log10(tab.turnover)
    tab['gap_to_next']=tab.log10_turnover.shift(-1)-tab.log10_turnover
    tab['ratio_to_next']=10**tab.gap_to_next
    n=len(tab); qmax=max(min_tail_size,int(math.ceil(max_fast_fraction*n)))
    aux=set(auxiliary_pools); rows=[]
    for i in range(n-1):
        g=float(tab.loc[i,'gap_to_next'])
        prev=float(tab.loc[i-1,'gap_to_next']) if i>0 else -np.inf
        nxt=float(tab.loc[i+1,'gap_to_next']) if i<n-2 else -np.inf
        q=n-i-1
        if q<=qmax and g>=prev and g>=nxt:
            raw=tab.loc[i+1:,'metabolite'].tolist()
            rows.append({'cut_after':tab.loc[i,'metabolite'],'gap_log10':g,'ratio':10**g,
                         'raw_fast_size':q,'raw_fast_pools':raw,
                         'eligible_fast_size':sum(m not in aux for m in raw),
                         'eligible_fast_pools':[m for m in raw if m not in aux]})
    cand=pd.DataFrame(rows).sort_values('gap_log10',ascending=False,ignore_index=True)
    return tab,cand


def primary_candidate(candidates: pd.DataFrame) -> pd.Series:
    if candidates.empty: raise ValueError('no top-tail local-gap candidate found')
    return candidates.iloc[0]

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument('csv',help='CSV with columns metabolite,turnover')
    ap.add_argument('--aux',nargs='*',default=[])
    args=ap.parse_args()
    df=pd.read_csv(args.csv)
    tab,cand=propose_fast_partitions(dict(zip(df.metabolite,df.turnover)),auxiliary_pools=args.aux)
    print(cand.to_string(index=False))
