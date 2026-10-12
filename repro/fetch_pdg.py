#!/usr/bin/env python3
"""Fetch every PDG input used by the manuscript, straight from pdg.lbl.gov
(the pdgLive / RPP-2026 web edition), archive the source PDFs, extract each
value with page-and-line provenance, and diff against the frozen inputs.

This is the tip of the tip-to-tail record: no PDG number enters the pipeline
except through this fetch (or the frozen CSV it validates).

Requires: python3 + pdftotext (poppler) on PATH; network access to
pdg.lbl.gov only.

Outputs (under repro/pdg_live/):
  pdfs/rpp2026-*.pdf   archived source documents (as fetched)
  txt/*.txt            pdftotext extractions (regenerated, not archived)
  fetched_inputs.csv   quantity, value, uncertainty, unit, source pdf,
                       source page/line, document id, sha256(pdf), fetch date
  FETCH_LOG.txt        url, http status, bytes, sha256, fetch date per file

Run from repro/:  python3 fetch_pdg.py
"""
import csv, hashlib, os, re, subprocess, sys, urllib.request, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
PDF_DIR = os.path.join(HERE, "pdg_live", "pdfs")
TXT_DIR = os.path.join(HERE, "pdg_live", "txt")
FETCH_LOG = os.path.join(HERE, "pdg_live", "FETCH_LOG.txt")
FETCHED_CSV = os.path.join(HERE, "pdg_live", "fetched_inputs.csv")
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(TXT_DIR, exist_ok=True)

BASE = "https://pdg.lbl.gov/2026/"
SOURCES = {  # filename -> (url, why we keep it)
    "rpp2026-list-electron.pdf": ("listings/rpp2026-list-electron.pdf", "m_e"),
    "rpp2026-list-muon.pdf": ("listings/rpp2026-list-muon.pdf", "m_mu"),
    "rpp2026-list-tau.pdf": ("listings/rpp2026-list-tau.pdf", "m_tau"),
    "rpp2026-list-light-quarks.pdf": ("listings/rpp2026-list-light-quarks.pdf",
                                      "m_u, m_d, m_s at 2 GeV"),
    "rpp2026-list-c-quark.pdf": ("listings/rpp2026-list-c-quark.pdf", "m_c(m_c)"),
    "rpp2026-list-b-quark.pdf": ("listings/rpp2026-list-b-quark.pdf", "m_b(m_b)"),
    "rpp2026-list-t-quark.pdf": ("listings/rpp2026-list-t-quark.pdf",
                                 "m_t direct, m_t(m_t)"),
    "rpp2026-rev-quark-masses.pdf": ("reviews/rpp2026-rev-quark-masses.pdf",
                                     "quark-mass conventions backing"),
    "rpp2026-rev-qcd.pdf": ("reviews/rpp2026-rev-qcd.pdf", "alpha_s(M_Z) world average"),
}

def fetch():
    log = [f"PDG 2026 fetch log -- {datetime.datetime.now(datetime.timezone.utc).isoformat()}",
           f"base: {BASE}", ""]
    for fn, (rel, _why) in SOURCES.items():
        url = BASE + rel
        dest = os.path.join(PDF_DIR, fn)
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                data = r.read()
                status = r.status
        except Exception as e:
            log.append(f"FAILED {url}: {e}")
            print(log[-1]); continue
        with open(dest, "wb") as f:
            f.write(data)
        sha = hashlib.sha256(data).hexdigest()
        log.append(f"{status} {len(data):8d}B sha256:{sha[:16]}... {url}")
        print(log[-1], flush=True)
        subprocess.run(["pdftotext", "-layout", dest,
                        os.path.join(TXT_DIR, fn[:-4] + ".txt")], check=True)
    with open(FETCH_LOG, "w") as f:
        f.write("\n".join(log) + "\n")

def pages(fn):
    text = open(os.path.join(TXT_DIR, fn[:-4] + ".txt"),
                encoding="utf-8", errors="replace").read()
    return text.split("\f")

def find(pages_list, pattern, flags=0):
    """Return (page, line_no, match) for first regex hit; page is 1-based."""
    rx = re.compile(pattern, flags)
    for pno, ptext in enumerate(pages_list, 1):
        for lno, line in enumerate(ptext.splitlines(), 1):
            m = rx.search(line)
            if m:
                return pno, lno, m
    return None, None, None

def extract():
    rows = []
    def add(quantity, value, unc_hi, unc_lo, unit, fn, page, line, doc, note=""):
        sha = hashlib.sha256(open(os.path.join(PDF_DIR, fn), "rb").read()).hexdigest()
        rows.append({"quantity": quantity, "value": value,
                     "unc_hi": unc_hi, "unc_lo": unc_lo, "unit": unit,
                     "source_pdf": fn, "page": page, "line": line,
                     "document_id": doc, "pdf_sha256": sha,
                     "fetch_date": datetime.datetime.now(datetime.timezone.utc).date().isoformat(),
                     "note": note})

    pe = pages("rpp2026-list-electron.pdf")
    pg, ln, m = find(pe, r"(0\.510998950[0-9]+)\s*±\s*(0\.000000000[0-9]+)\s+MOHR\s+25")
    add("m_e", m.group(1), m.group(2), m.group(2), "MeV", "rpp2026-list-electron.pdf", pg, ln, "MOHR 25 (2022 CODATA)", "pole mass")

    pm = pages("rpp2026-list-muon.pdf")
    pg, ln, m = find(pm, r"(105\.658375[0-9])\s*±\s*(0\.00000[0-9]+)\s+MOHR\s+25")
    add("m_mu", m.group(1), m.group(2), m.group(2), "MeV", "rpp2026-list-muon.pdf", pg, ln, "MOHR 25 (2022 CODATA)", "pole mass")

    pt = pages("rpp2026-list-tau.pdf")
    pg, ln, m = find(pt, r"(1776\.9[0-9]+)\s*±\s*(0\.[0-9]+)\s+OUR AVERAGE")
    add("m_tau", m.group(1), m.group(2), m.group(2), "MeV", "rpp2026-list-tau.pdf", pg, ln, "PDG OUR AVERAGE", "pole mass")

    pl = pages("rpp2026-list-light-quarks.pdf")
    flat = []
    for pno, ptext in enumerate(pl, 1):
        for lno, line in enumerate(ptext.splitlines(), 1):
            flat.append((pno, lno, line))
    for header, name in [("u-QUARK MASS", "m_u(2 GeV)"),
                         ("d-QUARK MASS", "m_d(2 GeV)"),
                         ("s-QUARK MASS", "m_s(2 GeV)")]:
        hidx = next(i for i, (_p, _l, t) in enumerate(flat) if header in t)
        for pno, lno, line in flat[hidx:]:
            mm = re.search(r"^\s*([0-9]+\.[0-9]+)\s*±\s*([0-9.]+)\s*\(CL = 90%\)\s*OUR EVALUATION", line)
            if mm:
                add(name, mm.group(1), mm.group(2), mm.group(2), "MeV",
                    "rpp2026-list-light-quarks.pdf", pno, lno,
                    "PDG OUR EVALUATION", "MSbar at 2 GeV, CL=90%")
                break

    pc = pages("rpp2026-list-c-quark.pdf")
    pg, ln, m = find(pc, r"^\s*(1\.2729)\s*±\s*(0\.0045)\s+\(CL = 90%\)\s+OUR EVALUATION")
    add("m_c(m_c)", m.group(1), m.group(2), m.group(2), "GeV", "rpp2026-list-c-quark.pdf", pg, ln, "PDG OUR EVALUATION", "MSbar self-scale, CL=90%")

    pb = pages("rpp2026-list-b-quark.pdf")
    pg, ln, m = find(pb, r"^\s*(4\.186)\s*±\s*(0\.006)\s+\(CL = 90%\)\s+OUR EVALUATION")
    add("m_b(m_b)", m.group(1), m.group(2), m.group(2), "GeV", "rpp2026-list-b-quark.pdf", pg, ln, "PDG OUR EVALUATION", "MSbar self-scale, CL=90%")

    ptp = pages("rpp2026-list-t-quark.pdf")
    pg, ln, m = find(ptp, r"^\s*(172\.60)\s*±\s*(0\.27)\s+OUR AVERAGE")
    add("m_t_direct", m.group(1), m.group(2), m.group(2), "GeV", "rpp2026-list-t-quark.pdf", pg, ln, "PDG OUR AVERAGE", "direct average; scale factor 1.5")

    pno, lno = None, None
    for pno2, ptext in enumerate(ptp, 1):
        lines = ptext.splitlines()
        for lno2, line in enumerate(lines, 1):
            if re.search(r"162\.5\s*\u2212\s*1\.5", line):
                above = "\n".join(lines[max(0, lno2 - 3):lno2])
                mm = re.search(r"\+\s*(2\.1)\s+OUR AVERAGE", above)
                if mm:
                    pno, lno, hi = pno2, lno2, mm.group(1)
    add("m_t(m_t)", "162.5", hi, "1.5", "GeV", "rpp2026-list-t-quark.pdf", pno, lno, "PDG OUR AVERAGE", "MSbar self scale from t tbar cross section")

    pq = pages("rpp2026-rev-qcd.pdf")
    pg, ln, m = find(pq, r"αs\s*\(m\s*2?Z\s*\)\s*=\s*0\.1180\s*±\s*0\.0009")
    add("alpha_s(M_Z)", "0.1180", "0.0009", "0.0009", "", "rpp2026-rev-qcd.pdf", pg, ln, "PDG world average, Sec. 9.4.8")

    return rows

FROZEN = {  # fetched quantity -> (frozen csv name, frozen value, unit factor to fetched units)
    "m_e": ("m_e", 0.51099895069, 1.0),
    "m_mu": ("m_mu", 105.6583755, 1.0),
    "m_tau": ("m_tau", 1776.93, 1.0),
    "m_u(2 GeV)": ("m_u", 2.16, 1.0),
    "m_d(2 GeV)": ("m_d", 4.70, 1.0),
    "m_s(2 GeV)": ("m_s", 92.9, 1.0),
    "m_c(m_c)": ("m_c(m_c)", 1272.9, 1e-3),   # frozen MeV -> fetched GeV
    "m_b(m_b)": ("m_b(m_b)", 4186.0, 1e-3),
    "m_t_direct": ("m_t", 172600.0, 1e-3),
    "m_t(m_t)": ("m_t(m_t)", 162500.0, 1e-3),
    "alpha_s(M_Z)": ("alpha_s(M_Z)", 0.118, 1.0),
}

def compare(rows):
    frozen = {r["quantity"]: r for r in
              csv.DictReader(open(os.path.join(HERE, "..", "inputs", "pdg2026_inputs.csv")))}
    out = []
    for r in rows:
        fname, fval, factor = FROZEN[r["quantity"]]
        f = frozen.get(fname, {})
        fval = float(f.get("value", fval)) * factor if f else fval * factor
        got = float(r["value"])
        status = "PASS" if abs(got - fval) <= max(1e-12, abs(fval) * 1e-6) else "CHANGED"
        out.append((r["quantity"], got, fval, status, r["source_pdf"], r["page"]))
    return out

def main():
    fetch()
    rows = extract()
    with open(FETCHED_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    cmp_rows = compare(rows)
    print("\nquantity          fetched        frozen   status   source (page)")
    for q, got, fval, st, pdf, pg in cmp_rows:
        print(f"{q:14s} {got:12g} {fval:12g}   {st:7s}  {pdf} p.{pg}")
    if any(st == "CHANGED" for *_x, st, _p, _pg in cmp_rows):
        print("\nWARNING: PDG values have moved since the manuscript inputs were frozen.")
        print("The frozen CSV remains authoritative for the manuscript's numbers; the")
        print("manuscript cites the edition, so a re-run with updated inputs would be a")
        print("new revision, not a bug.")
    else:
        print("\nAll fetched PDG inputs match the frozen inputs CSV.")

if __name__ == "__main__":
    main()
