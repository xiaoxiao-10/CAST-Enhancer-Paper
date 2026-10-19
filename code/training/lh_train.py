#!/usr/bin/env python
"""Locus-held-out retraining for the v2 audit (2026-10-17).

Holds out 20% of loci (stratified by chromosome and GC), retrains the three
model families on the remaining 80%, and evaluates on held-out cells restricted
to held-out loci (unseen) versus retained loci (seen).

Also trains a randomly initialised Nucleotide Transformer with identical
architecture and hyperparameters (the pretraining control).
"""
import argparse
import json
import os
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.allow_tf32 = True

from peft import LoraConfig, get_peft_model
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from transformers import AutoConfig, AutoModelForMaskedLM

BASE = "/home/vipuser/cast-v2-v0"
MODEL = ("/home/vipuser/CAST-Enhancer/data/raw/model_cache/"
         "nucleotide-transformer-v2-50m-multi-species-81b29e5")
USAGE = ("/home/vipuser/CAST-Enhancer/data/raw/fantom5/"
         "F5.hg38.enhancers.expression.usage.matrix.gz")
CELLS = {"K562": 0, "HepG2": 1, "GM12878": 2}
UNK = 3
VALID = {"chr18", "chr19", "chr14"}
TEST = {"chr3", "chr5", "chr11"}
SEEDS = [13, 37, 71, 101, 137]
DV = torch.device("cuda" if torch.cuda.is_available() else "cpu")
torch.set_num_threads(8)


def load_inputs():
    idx = pd.read_parquet(os.path.join(BASE, "v2_loci_index.parquet"))
    codes = np.load(os.path.join(BASE, "v2_seq_codes.npy"), mmap_mode="r")
    ids = np.load(os.path.join(BASE, "nt_ids_all.npy"), mmap_mode="r")
    am = np.load(os.path.join(BASE, "nt_am_all.npy"), mmap_mode="r")
    return idx, codes, ids, am


def build_dev_table(idx):
    import gzip
    with gzip.open(USAGE, "rt") as fh:
        header = [c.strip() for c in fh.readline().rstrip("\n").split("\t")]
    posmap = {c: i + 1 for i, c in enumerate(header)}
    devlibs = {"K562": ["CNhs12334", "CNhs12335", "CNhs12336"],
               "HepG2": ["CNhs12328", "CNhs12329", "CNhs12330"],
               "GM12878": ["CNhs12331", "CNhs12332", "CNhs12333"]}

    def build_lab(libs, rule):
        d = pd.read_csv(USAGE, sep="\t", compression="gzip",
                        usecols=[0] + [posmap[x] for x in libs], index_col=0)
        s = d.apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy().sum(axis=1)
        return pd.Series(np.where(s >= rule, 1, np.where(s == 0, 0, -1)), index=d.index)

    eid = idx["enhancer_id"].values
    chrom = dict(zip(idx["enhancer_id"], idx["chrom"]))
    eidx = {e: i for i, e in enumerate(eid)}
    rows = []
    for cell, libs in devlibs.items():
        lab = build_lab(libs, 2)
        for e in eid:
            v = lab.get(e, -1)
            if v < 0:
                continue
            ch = chrom[e]
            sp = "validation" if ch in VALID else ("test" if ch in TEST else "train")
            rows.append((eidx[e], CELLS[cell], int(v), sp))
    return pd.DataFrame(rows, columns=["r", "cidx", "y", "split"])


def build_holdout(idx, every=5):
    hol = np.zeros(len(idx), dtype=bool)
    for _, g in idx.groupby("chrom"):
        g = g.sort_values("gc")
        rows = g["row"].values.astype(int)
        take = rows[np.arange(len(rows)) % every == 0]
        hol[take] = True
    return hol


def comp_features(codes, n):
    out = np.zeros((n, 20), dtype=np.float32)
    for i in range(n):
        s = np.asarray(codes[i])
        s = s[s < 4]
        mono = np.bincount(s, minlength=4).astype(np.float32)
        mono /= max(mono.sum(), 1)
        dn = np.bincount(s[:-1] * 4 + s[1:], minlength=16).astype(np.float32)
        dn /= max(dn.sum(), 1)
        out[i, :4] = mono
        out[i, 4:] = dn
    return out


OH = np.zeros((5, 4), dtype=np.float32)
OH[:4] = np.eye(4, dtype=np.float32)
OH[4] = 0.25


class CNN(nn.Module):
    def __init__(s, cdim=32):
        super().__init__()
        s.emb = nn.Embedding(4, cdim)
        s.c1 = nn.Conv1d(4, 64, 11, padding=5)
        s.c2 = nn.Conv1d(64, 128, 7, padding=3)
        s.c3 = nn.Conv1d(128, 256, 5, padding=2)
        s.film = nn.Linear(cdim, 512)
        s.head = nn.Sequential(nn.Linear(512, 128), nn.GELU(),
                               nn.Dropout(0.2), nn.Linear(128, 1))
        s.p1 = nn.MaxPool1d(4)
        s.p2 = nn.MaxPool1d(4)
        s.p3 = nn.MaxPool1d(5)

    def forward(s, x, ci):
        h = torch.relu(s.c1(x)); h = s.p1(h)
        h = torch.relu(s.c2(h)); h = s.p2(h); h = s.c3(h)
        g, b = s.film(s.emb(ci)).chunk(2, dim=1)
        h = torch.relu(h * (1 + g.unsqueeze(-1)) + b.unsqueeze(-1))
        h = s.p3(h)
        h = torch.cat([h.max(dim=2).values, h.mean(dim=2)], dim=1)
        return s.head(h).squeeze(-1)


def to_tensor(codes, rows):
    x = OH[np.asarray(codes[rows])]
    return torch.from_numpy(np.ascontiguousarray(x.transpose(0, 2, 1))).float()


@torch.no_grad()
def cnn_score(net, codes, n, bs=256):
    net.eval()
    out = []
    for i in range(0, n, bs):
        k = min(bs, n - i)
        xb = to_tensor(codes, np.arange(i, i + k)).to(DV)
        cb = torch.full((k,), UNK, dtype=torch.long, device=DV)
        out.append(torch.sigmoid(net(xb, cb)).cpu().numpy())
    return np.concatenate(out)


class NTNet(nn.Module):
    def __init__(s, base, cdim=32):
        super().__init__()
        s.base = base
        s.emb = nn.Embedding(4, cdim)
        s.head = nn.Sequential(nn.Linear(base.config.hidden_size + cdim, 256),
                               nn.GELU(), nn.Dropout(0.2), nn.Linear(256, 1))

    def forward(s, ii, am, ci):
        o = s.base(input_ids=ii, attention_mask=am).last_hidden_state
        mk = am.unsqueeze(-1).float()
        pooled = (o * mk).sum(1) / mk.sum(1).clamp(min=1)
        return s.head(torch.cat([pooled, s.emb(ci)], dim=1)).squeeze(-1)


def nt_fresh(pretrained=True):
    if pretrained:
        base = AutoModelForMaskedLM.from_pretrained(MODEL, trust_remote_code=True).esm
    else:
        cfg = AutoConfig.from_pretrained(MODEL, trust_remote_code=True)
        base = AutoModelForMaskedLM.from_config(cfg, trust_remote_code=True).esm
    lcfg = LoraConfig(r=8, lora_alpha=16, lora_dropout=0.05,
                      target_modules=["query", "value"], task_type=None)
    return NTNet(get_peft_model(base, lcfg)).to(DV)


@torch.no_grad()
def nt_score(net, ids, am, n, bs=64):
    net.eval()
    out = []
    for i in range(0, n, bs):
        k = min(bs, n - i)
        ii = torch.from_numpy(ids[i:i + k]).long().to(DV)
        aa = torch.from_numpy(am[i:i + k]).long().to(DV)
        cc = torch.full((k,), UNK, dtype=torch.long, device=DV)
        out.append(torch.sigmoid(net(ii, aa, cc)).cpu().numpy())
    return np.concatenate(out)


def nt_train(seed, D, ids, am, epochs=4, bs=32, acc=2):
    torch.manual_seed(seed); np.random.seed(seed)
    net = nt_fresh(True)
    opt = torch.optim.AdamW([
        {"params": [p for n, p in net.named_parameters() if "lora_" in n], "lr": 2e-5},
        {"params": [p for n, p in net.named_parameters() if "lora_" not in n], "lr": 2e-4}],
        weight_decay=0.01)
    tr = D[D.split == "train"]; va = D[D.split == "validation"]
    pos = tr.y.sum()
    pw = torch.tensor([(len(tr) - pos) / max(pos, 1)], dtype=torch.float32, device=DV)
    lf = nn.BCEWithLogitsLoss(pos_weight=pw)
    R = tr.r.values.astype(np.int64); C = tr.cidx.values.astype(np.int64)
    Y = tr.y.values.astype(np.float32)
    best = -1.0; best_state = None
    for ep in range(1, epochs + 1):
        net.train(); perm = np.random.permutation(len(tr)); tot = 0.0; opt.zero_grad()
        for i in range(0, len(perm), bs):
            j = perm[i:i + bs]
            ii = torch.from_numpy(ids[R[j]]).long().to(DV)
            aa = torch.from_numpy(am[R[j]]).long().to(DV)
            cc = torch.from_numpy(C[j].copy()).long()
            m = torch.rand(len(cc)) < 0.3; cc[m] = UNK
            yy = torch.from_numpy(Y[j]).to(DV)
            loss = lf(net(ii, aa, cc.to(DV)), yy) / acc
            loss.backward(); tot += loss.item() * len(j) * acc
            if ((i // bs) + 1) % acc == 0:
                opt.step(); opt.zero_grad()
            del ii, aa, cc, yy
        opt.step(); opt.zero_grad()
        net.eval(); pv = []
        with torch.no_grad():
            Rv = va.r.values.astype(np.int64); Cv = va.cidx.values.astype(np.int64)
            for i in range(0, len(Rv), 64):
                sl = slice(i, i + 64)
                pv.append(torch.sigmoid(net(
                    torch.from_numpy(ids[Rv[sl]]).long().to(DV),
                    torch.from_numpy(am[Rv[sl]]).long().to(DV),
                    torch.from_numpy(Cv[sl]).long().to(DV))).cpu().numpy())
        ap = average_precision_score(va.y.values, np.concatenate(pv))
        print("  seed=%d ep=%d loss=%.4f val_auprc=%.4f" % (seed, ep, tot / len(perm), ap), flush=True)
        if ap > best:
            best = ap
            best_state = {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}
    net.load_state_dict(best_state)
    return net


def heldout_sets(idx):
    import gzip
    import json
    import re
    with gzip.open(USAGE, "rt") as fh:
        header = [c.strip() for c in fh.readline().rstrip("\n").split("\t")]
    posmap = {c: i + 1 for i, c in enumerate(header)}
    h = json.load(open(os.path.join(BASE, "heldout_cells_v2.json")))
    m = pd.read_csv(os.path.join(BASE, "all_library_map.csv"))
    m["library"] = m["library"].astype(str)
    m["celltype"] = m["sample"].astype(str).map(
        lambda s: re.sub(r",\s*donor\d+.*$", "", s).strip())
    primary = m[m["source"].str.contains("primary_cell")]
    eid = idx["enhancer_id"].values
    out = []
    for rec in h["selected"]:
        ct = rec["celltype"]
        libs = sorted(set(primary[primary.celltype == ct]["library"]) & set(header))
        need = int(np.ceil(2 / 3 * len(libs)))
        d = pd.read_csv(USAGE, sep="\t", compression="gzip",
                        usecols=[0] + [posmap[x] for x in libs], index_col=0)
        s = d.apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy().sum(axis=1)
        lab = pd.Series(np.where(s >= need, 1, np.where(s == 0, 0, -1)), index=d.index)
        rows = np.array([i for i, e in enumerate(eid) if lab.get(e, -1) >= 0], dtype=np.int64)
        y = np.array([lab[eid[i]] for i in rows], dtype=np.int64)
        out.append((ct, rows, y))
    return out


def evaluate(scores, ho_sets, hol, tag, rows_out):
    for ct, rows, y in ho_sets:
        for subset, mask in (("unseen", hol[rows]), ("seen", ~hol[rows])):
            r = rows[mask]; yy = y[mask]
            if yy.sum() < 10:
                continue
            p = scores[r]
            ap = average_precision_score(yy, p)
            rows_out.append(dict(model=tag, cell=ct, subset=subset, n=len(yy),
                                 pos=int(yy.sum()), prevalence=float(yy.mean()),
                                 auprc=float(ap), lift=float(ap / yy.mean()),
                                 auroc=float(roc_auc_score(yy, p))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/home/vipuser/cast-v2-locus-holdout")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--models", default="comp,cnn,nt,ntinit")
    args = ap.parse_args()
    models = set(args.models.split(","))
    os.makedirs(args.out, exist_ok=True)
    os.chdir(args.out)

    idx, codes, ids, am = load_inputs()
    n = len(idx)
    hol = build_holdout(idx)
    print("loci=%d holdout=%d (%.1f%%)" % (n, hol.sum(), 100 * hol.mean()), flush=True)
    np.save("LH_locus_holdout_mask.npy", hol)

    D = build_dev_table(idx)
    D = D[~hol[D.r.values]].reset_index(drop=True)
    if args.smoke:
        D = pd.concat([
            D[D.split == "train"].sample(min(3000, (D.split == "train").sum()), random_state=0),
            D[D.split == "validation"].sample(min(1000, (D.split == "validation").sum()), random_state=0),
            D[D.split == "test"].sample(min(1000, (D.split == "test").sum()), random_state=0),
        ]).reset_index(drop=True)
    print("dev rows after holdout:", len(D), D.groupby("split").size().to_dict(), flush=True)
    pd.DataFrame([dict(split=k, n=v) for k, v in D.groupby("split").size().items()]).to_csv(
        "LH_dev_split_counts.tsv", sep="\t", index=False)

    ho_sets = heldout_sets(idx)
    rows_out = []
    seeds = [13] if args.smoke else SEEDS

    # composition
    if "comp" in models:
        t0 = time.time()
        X = comp_features(codes, n)
        tr = D[D.split == "train"]
        clf = LogisticRegression(max_iter=3000, C=1.0, class_weight="balanced").fit(
            X[tr.r.values], tr.y.values)
        p_comp = clf.predict_proba(X)[:, 1]
        np.save("LH_scores_comp.npy", p_comp)
        print("composition done %.1fs" % (time.time() - t0), flush=True)

    # CNN
    cnn_scores = []
    for seed in (seeds if "cnn" in models else []):
        torch.manual_seed(seed); np.random.seed(seed)
        net = CNN().to(DV)
        tr = D[D.split == "train"]; va = D[D.split == "validation"]
        pos = tr.y.sum()
        pw = torch.tensor([(len(tr) - pos) / max(pos, 1)], dtype=torch.float32, device=DV)
        lf = nn.BCEWithLogitsLoss(pos_weight=pw)
        opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=0.01)
        R = tr.r.values.astype(np.int64); C = tr.cidx.values.astype(np.int64)
        Y = tr.y.values.astype(np.float32)
        best = -1; best_state = None; stale = 0
        epochs = 1 if args.smoke else 15
        for ep in range(1, epochs + 1):
            net.train(); perm = np.random.permutation(len(tr)); tot = 0.0
            for i in range(0, len(perm), 64):
                j = perm[i:i + 64]
                xb = to_tensor(codes, R[j]).to(DV)
                cb = torch.from_numpy(C[j].copy()).long()
                m = torch.rand(len(cb)) < 0.3; cb[m] = UNK
                yb = torch.from_numpy(Y[j]).to(DV)
                opt.zero_grad(); loss = lf(net(xb, cb.to(DV)), yb)
                loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 1.0); opt.step()
                tot += loss.item() * len(j); del xb, cb, yb
            va_r = va.r.values.astype(np.int64); va_c = va.cidx.values.astype(np.int64)
            pv = []
            net.eval()
            with torch.no_grad():
                for i in range(0, len(va_r), 256):
                    xb = to_tensor(codes, va_r[i:i + 256]).to(DV)
                    cb = torch.from_numpy(va_c[i:i + 256]).long().to(DV)
                    pv.append(torch.sigmoid(net(xb, cb)).cpu().numpy())
            val_ap = average_precision_score(va.y.values, np.concatenate(pv))
            print("  CNN seed=%d ep=%d loss=%.4f val=%.4f" % (seed, ep, tot / len(perm), val_ap), flush=True)
            if val_ap > best + 1e-4:
                best = val_ap; best_state = {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}; stale = 0
            else:
                stale += 1
                if stale >= 3: break
        net.load_state_dict(best_state); torch.save(best_state, "LH_cnn_seed_%d.pt" % seed)
        s = cnn_score(net, codes, n if not args.smoke else min(n, 2000))
        np.save("LH_scores_cnn_seed_%d.npy" % seed, s)
        cnn_scores.append(s); del net; torch.cuda.empty_cache()
    p_cnn = np.mean(cnn_scores, axis=0) if cnn_scores else None
    if not args.smoke and p_cnn is not None:
        np.save("LH_scores_cnn_ensemble.npy", p_cnn)

    # NT pretrained + random init
    for tag, pretrained in (("nt", True), ("ntinit", False)):
        if tag not in models:
            continue
        all_scores = []
        for seed in seeds:
            torch.manual_seed(seed); np.random.seed(seed)
            if pretrained:
                net = nt_train(seed, D, ids, am, epochs=1 if args.smoke else 4)
            else:
                net = nt_fresh(False)
                opt = torch.optim.AdamW([
                    {"params": [p for nm, p in net.named_parameters() if "lora_" in nm], "lr": 2e-5},
                    {"params": [p for nm, p in net.named_parameters() if "lora_" not in nm], "lr": 2e-4}],
                    weight_decay=0.01)
                tr = D[D.split == "train"]; va = D[D.split == "validation"]
                pos = tr.y.sum()
                pw = torch.tensor([(len(tr) - pos) / max(pos, 1)], dtype=torch.float32, device=DV)
                lf = nn.BCEWithLogitsLoss(pos_weight=pw)
                R = tr.r.values.astype(np.int64); C = tr.cidx.values.astype(np.int64)
                Y = tr.y.values.astype(np.float32)
                best = -1; best_state = None
                for ep in range(1, (1 if args.smoke else 4) + 1):
                    net.train(); perm = np.random.permutation(len(tr)); tot = 0.0; opt.zero_grad()
                    for i in range(0, len(perm), 32):
                        j = perm[i:i + 32]
                        ii = torch.from_numpy(ids[R[j]]).long().to(DV)
                        aa = torch.from_numpy(am[R[j]]).long().to(DV)
                        cc = torch.from_numpy(C[j].copy()).long()
                        m = torch.rand(len(cc)) < 0.3; cc[m] = UNK
                        yy = torch.from_numpy(Y[j]).to(DV)
                        loss = lf(net(ii, aa, cc.to(DV)), yy) / 2
                        loss.backward(); tot += loss.item() * len(j) * 2
                        if ((i // 32) + 1) % 2 == 0: opt.step(); opt.zero_grad()
                        del ii, aa, cc, yy
                    opt.step(); opt.zero_grad()
                    net.eval(); pv = []
                    with torch.no_grad():
                        Rv = va.r.values.astype(np.int64); Cv = va.cidx.values.astype(np.int64)
                        for i in range(0, len(Rv), 64):
                            sl = slice(i, i + 64)
                            pv.append(torch.sigmoid(net(
                                torch.from_numpy(ids[Rv[sl]]).long().to(DV),
                                torch.from_numpy(am[Rv[sl]]).long().to(DV),
                                torch.from_numpy(Cv[sl]).long().to(DV))).cpu().numpy())
                    val_ap = average_precision_score(va.y.values, np.concatenate(pv))
                    print("  %s seed=%d ep=%d loss=%.4f val=%.4f" % (tag, seed, ep, tot / len(perm), val_ap), flush=True)
                    if val_ap > best:
                        best = val_ap; best_state = {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}
                net.load_state_dict(best_state)
            torch.save(net.state_dict(), "LH_%s_seed_%d.pt" % (tag, seed))
            s = nt_score(net, ids, am, n if not args.smoke else min(n, 2000))
            np.save("LH_scores_%s_seed_%d.npy" % (tag, seed), s)
            all_scores.append(s); del net; torch.cuda.empty_cache()
        if not args.smoke:
            np.save("LH_scores_%s_ensemble.npy" % tag, np.mean(all_scores, axis=0))

    if args.smoke:
        print("SMOKE OK: comp, CNN, NT pretrained and random-init all trained/scored", flush=True)
        return

    for tag, fname in (("comp", "LH_scores_comp.npy"),
                       ("cnn", "LH_scores_cnn_ensemble.npy"),
                       ("nt", "LH_scores_nt_ensemble.npy"),
                       ("ntinit", "LH_scores_ntinit_ensemble.npy")):
        if os.path.exists(fname):
            evaluate(np.load(fname), ho_sets, hol, tag, rows_out)
    R = pd.DataFrame(rows_out)
    R.to_csv("LH_results.tsv", sep="\t", index=False)
    print(R.groupby(["model", "subset"])[["auprc", "lift"]].mean().to_string(), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
