"""Run all eight sample exports and descriptive two-sample age-group summaries.

Usage: python scripts/run_liana.py EXPORT_ROOT NEW_OUTPUT_DIRECTORY config/samples.tsv
Original defaults are retained. Missing distal estimates are filled with zero, as in
source notebooks; this convention does not establish that signaling is absent.
"""
from pathlib import Path
import argparse
import importlib.metadata
import json
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import liana as li


def load_export(folder):
    folder = Path(folder)
    x = sio.mmread(folder / "X.mtx").tocsr().T.tocsr()
    genes = pd.read_csv(folder / "genes.txt", header=None)[0].astype(str).tolist()
    cells = pd.read_csv(folder / "cells.txt", header=None)[0].astype(str).tolist()
    obs = pd.read_csv(folder / "obs.csv").set_index("cell_id")
    if len(set(cells)) != len(cells) or len(set(genes)) != len(genes):
        raise ValueError("Duplicated cell or gene identifiers")
    if not set(cells).issubset(obs.index):
        raise ValueError("Expression cells absent from metadata")
    if x.shape != (len(cells), len(genes)):
        raise ValueError("Matrix and identifier lengths differ")
    a = sc.AnnData(x, obs=obs.loc[cells].copy())
    a.var_names, a.obs_names = genes, cells
    if "Core" not in set(a.obs.CCI_group):
        raise ValueError("Core group is absent")
    return a


def neighbor_distal_keep(res):
    required = {"target", "source", "ligand_complex", "receptor_complex", "lrscore"}
    if not required.issubset(res.columns):
        raise ValueError(f"LIANA result lacks required source-workflow columns: {required-set(res.columns)}")
    df = res.copy()
    df["lr"] = df.ligand_complex.astype(str) + "→" + df.receptor_complex.astype(str)
    nb = df[(df.target == "Core") & df.source.str.startswith("Neighbor_")].copy()
    ds = df[(df.target == "Core") & df.source.str.startswith("Distal_")].copy()
    for frame, prefix in ((nb, "Neighbor_"), (ds, "Distal_")):
        frame["ct"] = frame.source.str.slice(len(prefix))
    nb = nb[~nb.ct.isin(["Erythroid cell", "Unknown"])]
    ds = ds[~ds.ct.isin(["Erythroid cell", "Unknown"])]
    nb = nb.groupby(["ct", "lr"], as_index=False).lrscore.mean().rename(columns={"lrscore": "mean_neighbor"})
    ds = ds.groupby(["ct", "lr"], as_index=False).lrscore.mean().rename(columns={"lrscore": "mean_distal"})
    m = nb.merge(ds, on=["ct", "lr"], how="left")
    m["distal_estimate_missing"] = m.mean_distal.isna()
    m["mean_distal"] = m.mean_distal.fillna(0.0)
    m["delta"] = m.mean_neighbor - m.mean_distal
    return m[m.delta > 0].copy()


def combine_two(k1, k2, prefix):
    cols = ["ct", "lr", "mean_neighbor"]
    x = k1[cols].rename(columns={"mean_neighbor": f"{prefix}_1"})
    y = k2[cols].rename(columns={"mean_neighbor": f"{prefix}_2"})
    z = x.merge(y, on=["ct", "lr"], how="outer")
    c = [f"{prefix}_1", f"{prefix}_2"]
    z[c] = z[c].fillna(0.0)
    z[f"{prefix}_present_n"] = (z[c] > 0).sum(axis=1)
    z[f"{prefix}_mean"] = z[c].mean(axis=1)
    return z


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("export_root", type=Path); p.add_argument("output", type=Path)
    p.add_argument("manifest", type=Path)
    args = p.parse_args()
    if args.output.exists():
        raise SystemExit("Choose a new output directory")
    m = pd.read_csv(args.manifest, sep="\t")
    if len(m) != 8 or m["sample"].duplicated().any():
        raise ValueError("Expected eight unique samples")
    for tag in m.liana_tag:
        for name in ("X.mtx", "genes.txt", "cells.txt", "obs.csv"):
            if not (args.export_root / tag / name).is_file():
                raise FileNotFoundError(args.export_root / tag / name)
    args.output.mkdir(parents=True)
    versions = {x: importlib.metadata.version(x) for x in ["liana", "scanpy", "anndata", "numpy", "pandas", "scipy"]}
    (args.output/"python_versions.json").write_text(json.dumps(versions,indent=2))
    keep = {}
    for row in m.itertuples(index=False):
        a = load_export(args.export_root / row.liana_tag)
        li.mt.rank_aggregate(a, groupby="CCI_group", resource_name="mouseconsensus", use_raw=False, verbose=True)
        res = a.uns["liana_res"].copy()
        res.to_csv(args.output/f"{row.sample}_liana_full.csv",index=False)
        k = neighbor_distal_keep(res)
        k.to_csv(args.output/f"{row.sample}_neighbor_gt_distal.csv",index=False)
        keep[row.sample] = k
    for tissue in ["Heart", "Adipose"]:
        a = combine_two(keep[f"{tissue}_Aged_1"],keep[f"{tissue}_Aged_2"],"Aged")
        y = combine_two(keep[f"{tissue}_Young_1"],keep[f"{tissue}_Young_2"],"Young")
        z = a.merge(y,on=["ct","lr"],how="outer")
        numeric_cols = [c for c in z if c not in ["ct","lr"]]
        z[numeric_cols] = z[numeric_cols].fillna(0.0)
        z = z[(z.Aged_present_n >= 2) | (z.Young_present_n >= 2)].copy()
        z["max_mean"] = z[["Aged_mean","Young_mean"]].max(axis=1)
        z.to_csv(args.output/f"{tissue}_age_group_summary.csv",index=False)

if __name__ == "__main__":
    main()
