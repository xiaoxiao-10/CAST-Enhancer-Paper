#!/usr/bin/env python
"""Build the Word draft from paper/manuscript_v3.md."""
import os
import re

import pandas as pd
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Inches, Pt, RGBColor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MD = os.path.join(ROOT, "manuscript", "en", "manuscript_v3.md")
FIG = os.path.join(ROOT, "figures", "main")
TAB = os.path.join(ROOT, "tables")
LH = os.path.join(ROOT, "results", "locus_holdout")
OUT = os.path.join(ROOT, "manuscript", "en", "CAST-Enhancer_v3_manuscript_draft.docx")

MODEL = {"comp": "Composition", "cnn": "DNA-CNN", "nt": "NT-LoRA"}
SUBSET = {"unseen": "Unseen", "seen": "Retained"}
DACT = {"all": "All", "dev_active": "Active", "dev_inactive": "Inactive"}

SHORT_LEGENDS = {
    "Fig. 1": "A genomic locus anchors the locus-held-out, cell-held-out design. (a) A representative FANTOM5 enhancer; ATAC tracks are shown for context only and are not model inputs. The locus is ATAC-open in all three development cells and CAGE-active only in K562. (b) Activity of 24 loci across 14 cell types, split into retained and held-out loci. (c) The three model families compared on identical inputs and splits. (d) The four evaluation strata.",
    "Fig. 2": "Stronger composition baselines remove the apparent NT advantage. Macro AUPRC for logistic composition, composition-XGBoost, motif-XGBoost, composition+motif-XGBoost and NT-LoRA across four locus strata. Δ is the macro NT − composition-XGBoost difference with its 95% paired block-bootstrap interval.",
    "Fig. 3": "The NT is indistinguishable from composition-XGBoost on unseen development-active loci. (a) Per-cell paired AUPRC; x is composition-XGBoost and y is NT-LoRA; the dashed diagonal marks equal performance. (b) Macro NT − composition-XGBoost ΔAUPRC with 95% paired block-bootstrap intervals across four strata; the dashed line marks zero.",
    "Fig. 4": "Pretrained initialization is required for the logistic-baseline contrast on unseen loci under this budget. (a) Macro AUPRC on retained and unseen development-active loci; thin lines are individual cell types. (b) Per-cell difference between pretrained and random-init NT on retained (x) and unseen (y) loci; all 11 cells lie above the diagonal. n = 11; exact intervals in Table S4.",
}

ALT_TEXT = {
    "fig1_design": "Design of the locus-held-out, cell-held-out experiment.",
    "fig2_four_strata": "Macro AUPRC for five composition, motif and NT models across four locus strata.",
    "fig3_per_cell": "Per-cell scatter and macro forest plot of NT-LoRA versus composition-XGBoost.",
    "fig4_control": "Pretraining versus random initialization on retained and unseen loci.",
}


def clean(text):
    text = re.sub(r"<!--.*?-->", "", text)
    text = text.replace("`", "")
    return re.sub(r"\s+", " ", text).strip()


def add_rich(p, text, size=None, color=None):
    parts = re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*)", clean(text))
    for part in parts:
        if not part:
            continue
        bold = part.startswith("**") and part.endswith("**")
        italic = part.startswith("*") and part.endswith("*") and not bold
        run = p.add_run(part[2:-2] if bold else part[1:-1] if italic else part)
        run.bold = bold
        run.italic = italic
        if size:
            run.font.size = Pt(size)
        if color:
            run.font.color.rgb = color
    return p


def set_cell_shading(cell, color):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), color)
    tcPr.append(shd)


def set_table_borders(table):
    """Apply a journal-style three-line table border."""
    tbl = table._tbl
    tblPr = tbl.tblPr
    for old in tblPr.findall(qn("w:tblBorders")):
        tblPr.remove(old)
    borders = OxmlElement("w:tblBorders")
    for edge, val, sz in (
        ("top", "single", "8"),
        ("bottom", "single", "8"),
        ("left", "nil", "0"),
        ("right", "nil", "0"),
        ("insideH", "nil", "0"),
        ("insideV", "nil", "0"),
    ):
        el = OxmlElement("w:" + edge)
        el.set(qn("w:val"), val)
        el.set(qn("w:sz"), sz)
        el.set(qn("w:space"), "0")
        el.set(qn("w:color"), "000000")
        borders.append(el)
    tblPr.append(borders)
    mar = OxmlElement("w:tblCellMar")
    for side in ("top", "left", "bottom", "right"):
        el = OxmlElement("w:" + side)
        el.set(qn("w:w"), "60")
        el.set(qn("w:type"), "dxa")
        mar.append(el)
    tblPr.append(mar)


def set_cell_edge(cell, edge, sz=6, color="000000"):
    """Set one visible edge on a cell, used for the header/footer rules."""
    tcPr = cell._tc.get_or_add_tcPr()
    tcBorders = tcPr.find(qn("w:tcBorders"))
    if tcBorders is None:
        tcBorders = OxmlElement("w:tcBorders")
        tcPr.append(tcBorders)
    old = tcBorders.find(qn("w:" + edge))
    if old is not None:
        tcBorders.remove(old)
    el = OxmlElement("w:" + edge)
    el.set(qn("w:val"), "single")
    el.set(qn("w:sz"), str(sz))
    el.set(qn("w:space"), "0")
    el.set(qn("w:color"), color)
    tcBorders.append(el)


def add_table(doc, caption, df, widths, font=8.2, header_fill="D9E2F3", page_break_before=False):
    cap = doc.add_paragraph()
    if page_break_before:
        cap.paragraph_format.page_break_before = True
    cap.paragraph_format.space_before = Pt(10)
    cap.paragraph_format.space_after = Pt(4)
    cap.paragraph_format.keep_with_next = True
    add_rich(cap, caption, size=9)
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    set_table_borders(t)
    hdr = t.rows[0]
    trPr = hdr._tr.get_or_add_trPr()
    tblHeader = OxmlElement("w:tblHeader")
    tblHeader.set(qn("w:val"), "true")
    trPr.append(tblHeader)
    for i, col in enumerate(df.columns):
        cell = hdr.cells[i]
        cell.width = Cm(widths[i])
        p = cell.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(2)
        p.paragraph_format.space_after = Pt(2)
        r = p.add_run(str(col))
        r.bold = True
        r.font.size = Pt(font)
        set_cell_shading(cell, header_fill)
    for ri, row in df.iterrows():
        cells = t.add_row().cells
        for i, val in enumerate(row):
            cell = cells[i]
            cell.width = Cm(widths[i])
            p = cell.paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT if i <= 1 else WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_before = Pt(1)
            p.paragraph_format.space_after = Pt(1)
            r = p.add_run(str(val))
            r.font.size = Pt(font)
            if ri % 2 == 1:
                set_cell_shading(cell, "F7F9FC")
    doc.add_paragraph().paragraph_format.space_after = Pt(2)
    return t


def add_page_number(section):
    p = section.footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    fld = OxmlElement("w:fldChar"); fld.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText"); instr.set(qn("xml:space"), "preserve"); instr.text = "PAGE"
    end = OxmlElement("w:fldChar"); end.set(qn("w:fldCharType"), "end")
    run._r.append(fld); run._r.append(instr); run._r.append(end)
    run.font.size = Pt(9)


def move_tables_before_discussion(doc):
    """Move the generated Tables block from the document end to before Discussion."""
    body = doc.element.body
    paragraphs = [p for p in doc.paragraphs if p.text.strip()]
    p_disc = next((p for p in paragraphs if p.text.strip().startswith("3. Discussion")), None)
    p_tables = next((p for p in paragraphs if p.text.strip() == "Tables"), None)
    if p_disc is None or p_tables is None:
        return

    children = list(body)
    table_elements = [el for el in children if el.tag == qn("w:tbl")]
    if len(table_elements) < 2:
        return

    start = children.index(p_tables._p)
    end = children.index(table_elements[-1])
    block = children[start:end + 1]
    for element in block:
        p_disc._p.addprevious(element)
    p_disc.paragraph_format.page_break_before = False


def parse_legends(text):
    m = re.search(r"## Figure legends(.*?)\n## ", text, re.S)
    section = m.group(1) if m else ""
    parts = re.split(r"\*\*(Fig\. \d+) \|", section)
    legends = {}
    for i in range(1, len(parts), 2):
        key = parts[i]
        chunk = parts[i + 1].split("**Table")[0]
        title, _, body = chunk.partition("**")
        legends[key] = clean(title + " " + body)
    return legends


def main():
    text = open(MD, encoding="utf-8").read()
    legends = parse_legends(text)
    doc = Document()
    sec = doc.sections[0]
    sec.page_width = Cm(21.0); sec.page_height = Cm(29.7)
    sec.left_margin = Cm(2.4); sec.right_margin = Cm(2.4)
    sec.top_margin = Cm(2.3); sec.bottom_margin = Cm(2.3)
    add_page_number(sec)

    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(12)
    normal.paragraph_format.line_spacing = 1.25
    normal.paragraph_format.space_after = Pt(6)

    # Remove the default bottom border from the Title style.
    title_style = doc.styles["Title"]
    t_pPr = title_style.element.get_or_add_pPr()
    for bdr in t_pPr.findall(qn("w:pBdr")):
        t_pPr.remove(bdr)

    for level, size, before, after, italic in ((1, 14, 14, 6, False), (2, 12.5, 10, 4, False), (3, 12, 8, 3, True)):
        st = doc.styles["Heading %d" % level]
        st.font.name = "Times New Roman"; st.font.size = Pt(size); st.font.bold = True
        st.font.italic = italic; st.font.color.rgb = RGBColor(0, 0, 0)
        st.paragraph_format.space_before = Pt(before)
        st.paragraph_format.space_after = Pt(after)
        st.paragraph_format.keep_with_next = True

    lines = text.splitlines()
    i = 0
    title_done = False
    skip = None
    table1_done = table2_done = False
    para = []
    bullet_para = None
    in_refs = False
    in_21 = False
    fig1_inserted = False

    def flush():
        nonlocal para
        if para:
            p = doc.add_paragraph()
            add_rich(p, " ".join(para))
            para = []

    def add_fig(stem):
        png = os.path.join(FIG, stem + ".png")
        key = {"fig1_design": "Fig. 1", "fig2_four_strata": "Fig. 2",
               "fig3_per_cell": "Fig. 3", "fig4_control": "Fig. 4"}.get(stem)
        width = {
            "fig1_design": 6.35,
            "fig2_four_strata": 6.35,
            "fig3_per_cell": 6.35,
            "fig4_sensitivity": 6.0,
            "fig4_control": 6.35,
        }.get(stem, 6.0)
        p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(6); p.paragraph_format.space_after = Pt(3)
        p.paragraph_format.keep_with_next = True
        p.add_run().add_picture(png, width=Inches(width))
        shape = doc.inline_shapes[-1]
        docPr = shape._inline.docPr
        docPr.set("descr", ALT_TEXT.get(stem, ""))
        docPr.set("title", ALT_TEXT.get(stem, ""))
        if key and key in SHORT_LEGENDS:
            cp = doc.add_paragraph()
            cp.paragraph_format.space_after = Pt(10)
            cp.paragraph_format.line_spacing = 1.1
            cp.paragraph_format.keep_together = True
            add_rich(cp, "**%s |** %s" % (key, SHORT_LEGENDS[key]), size=8.5)

    while i < len(lines):
        line = lines[i].rstrip()
        if line.startswith("# ") and not title_done:
            p = doc.add_paragraph(style="Title")
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p_pPr = p._p.get_or_add_pPr()
            for bdr in p_pPr.findall(qn("w:pBdr")):
                p_pPr.remove(bdr)
            r = p.add_run(line[2:].strip())
            r.font.size = Pt(16); r.bold = True; r.font.name = "Times New Roman"
            r.font.color.rgb = RGBColor(0, 0, 0)
            sub = doc.add_paragraph()
            sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
            rs = sub.add_run("Draft manuscript | 19 October 2026")
            rs.italic = True; rs.font.size = Pt(9); rs.font.name = "Times New Roman"
            title_done = True
            i += 1
            continue
        if line.startswith("## "):
            flush()
            bullet_para = None
            name = line[3:].strip()
            if name.startswith("Figure legends") or name.startswith("起草状态"):
                skip = name
                i += 1
                continue
            if name.startswith("References"):
                skip = None
                in_refs = True
                doc.add_heading("References", level=1)
                i += 1
                continue
            skip = None
            in_refs = False
            doc.add_heading(name, level=1)
            i += 1
            continue
        if skip:
            i += 1
            continue
        if line.startswith("### "):
            flush()
            bullet_para = None
            name = line[4:].strip()
            in_21 = name.startswith("2.1")
            doc.add_heading(name, level=2)
            i += 1
            continue
        if line.startswith("#### "):
            flush()
            bullet_para = None
            doc.add_heading(line[5:].strip(), level=3)
            i += 1
            continue
        if line.startswith("!["):
            flush()
            bullet_para = None
            m = re.search(r"\((figures_v3/([a-z0-9_]+)\.svg)\)", line)
            if m:
                if not (m.group(2) == "fig1_design" and fig1_inserted):
                    add_fig(m.group(2))
            i += 1
            continue
        if line.strip() == "":
            had = bool(para)
            flush()
            bullet_para = None
            if in_21 and had and not fig1_inserted:
                add_fig("fig1_design")
                fig1_inserted = True
                in_21 = False
            i += 1
            continue
        if line.strip() == "---":
            i += 1
            continue
        if in_refs and not re.match(r"^\d+\.\s+", line):
            i += 1
            continue
        if re.match(r"^\d+\.\s+", line):
            flush()
            bullet_para = None
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Cm(0.6)
            p.paragraph_format.first_line_indent = Cm(-0.6)
            p.paragraph_format.space_after = Pt(3)
            p.paragraph_format.line_spacing = 1.1
            add_rich(p, line, size=10)
            i += 1
            continue
        if line.startswith("- "):
            flush()
            bullet_para = doc.add_paragraph(style="List Bullet")
            bullet_para.paragraph_format.space_after = Pt(4)
            add_rich(bullet_para, line[2:].strip())
            i += 1
            continue
        if bullet_para is not None:
            bullet_para.add_run(" ")
            add_rich(bullet_para, line.strip())
            i += 1
            continue
        para.append(line.strip())
        i += 1
    flush()

    h = doc.add_heading("Tables", level=1)
    h.paragraph_format.page_break_before = False

    # ---------------- Table 1 ----------------
    cap = doc.add_paragraph()
    cap.paragraph_format.space_before = Pt(6)
    cap.paragraph_format.space_after = Pt(4)
    cap.paragraph_format.line_spacing = 1.15
    cap.paragraph_format.keep_with_next = True
    add_rich(cap, "**Table 1 | Cohort and label counts.** Development cells after the 20% locus holdout, and the 11 held-out primary cell types.", size=9)
    df = pd.read_csv(os.path.join(TAB, "table1_cohort.tsv"), sep="\t")
    widths1 = [3.4, 2.4, 1.8, 1.8, 1.8, 1.8, 2.2]
    headers1 = ["Cell type", "Split", "Libraries", "Threshold", "Positives", "Negatives", "Prevalence"]
    t = doc.add_table(rows=1, cols=len(headers1))
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    set_table_borders(t)
    hdr = t.rows[0]
    trPr = hdr._tr.get_or_add_trPr()
    th = OxmlElement("w:tblHeader"); th.set(qn("w:val"), "true"); trPr.append(th)
    for i, text in enumerate(headers1):
        cell = hdr.cells[i]; cell.width = Cm(widths1[i])
        p = cell.paragraphs[0]; p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(2); p.paragraph_format.space_after = Pt(2)
        p.paragraph_format.line_spacing = 1.0
        r = p.add_run(text); r.bold = True; r.font.size = Pt(7.6)
        r.font.color.rgb = RGBColor(0, 0, 0)
        set_cell_edge(cell, "top", 8)
        set_cell_edge(cell, "bottom", 6)

    def group_row(text):
        row = t.add_row()
        cell = row.cells[0].merge(row.cells[-1])
        p = cell.paragraphs[0]
        p.paragraph_format.space_before = Pt(7)
        p.paragraph_format.space_after = Pt(2)
        p.paragraph_format.line_spacing = 1.0
        r = p.add_run(text); r.bold = True; r.font.size = Pt(7.8)

    def data_row(vals, shade):
        cells = t.add_row().cells
        for i, (v, w) in enumerate(zip(vals, widths1)):
            cell = cells[i]; cell.width = Cm(w)
            p = cell.paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT if i == 0 else WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_before = Pt(1)
            p.paragraph_format.space_after = Pt(1)
            p.paragraph_format.line_spacing = 1.0
            r = p.add_run(str(v)); r.font.size = Pt(7.4)

    group_row("Development cells (after 20% locus holdout)")
    dev = df[df.block == "development"]
    for k, row in enumerate(dev.itertuples()):
        data_row([row.cell, row.split, row.n_libraries, row.threshold, row.positives,
                  row.negatives, "%.2f%%" % (100 * row.prevalence)], k % 2 == 1)
    group_row("Held-out primary cell types (labels read once)")
    ho = df[df.block == "held-out"]
    for k, row in enumerate(ho.itertuples()):
        data_row([row.cell, "held-out", row.n_libraries, row.threshold, row.positives,
                  row.negatives, "%.2f%%" % (100 * row.prevalence)], k % 2 == 1)
    for cell in t.rows[-1].cells:
        set_cell_edge(cell, "bottom", 8)
    doc.add_paragraph().paragraph_format.space_after = Pt(2)

    # ---------------- Table 2 ----------------
    cap = doc.add_paragraph()
    cap.paragraph_format.space_before = Pt(6)
    cap.paragraph_format.space_after = Pt(4)
    cap.paragraph_format.line_spacing = 1.15
    cap.paragraph_format.keep_with_next = True
    add_rich(cap, "**Table 2 | Four-stratum performance.** AUPRC with mean per-cell lift in parentheses; n and positives are per-cell means, and the unseen-inactive row uses the 9 cell types with at least 10 positives; Δ is the NT-LoRA minus logistic composition point estimate with 95% CI.", size=9)
    perf = pd.read_csv(os.path.join(TAB, "table2_four_strata.tsv"), sep="\t")
    macro = pd.read_csv(os.path.join(LH, "LH_macro_point.tsv"), sep="\t")
    order = [("unseen", "all"), ("unseen", "dev_active"),
             ("seen", "dev_active"), ("unseen", "dev_inactive")]
    cols2 = ["Locus set", "Development activity", "n / positives", "Composition (logistic)",
             "DNA-CNN", "NT-LoRA", "Δ AUPRC (95% CI)"]
    widths2 = [1.8, 2.3, 2.1, 2.2, 2.2, 2.2, 3.0]
    t2 = doc.add_table(rows=1, cols=len(cols2))
    t2.alignment = WD_TABLE_ALIGNMENT.CENTER
    t2.autofit = False
    set_table_borders(t2)
    hdr = t2.rows[0]
    trPr = hdr._tr.get_or_add_trPr()
    th = OxmlElement("w:tblHeader"); th.set(qn("w:val"), "true"); trPr.append(th)
    for i, text in enumerate(cols2):
        cell = hdr.cells[i]; cell.width = Cm(widths2[i])
        p = cell.paragraphs[0]; p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(3); p.paragraph_format.space_after = Pt(3)
        p.paragraph_format.line_spacing = 1.0
        r = p.add_run(text); r.bold = True; r.font.size = Pt(7.6)
        r.font.color.rgb = RGBColor(0, 0, 0)
        set_cell_edge(cell, "top", 8)
        set_cell_edge(cell, "bottom", 6)
    for k, (sub, dact) in enumerate(order):
        cells = t2.add_row().cells
        base_row = perf[(perf.subset == sub) & (perf.dact == dact)].iloc[0]
        vals = [SUBSET[sub], DACT[dact],
                "%s / %d" % (format(int(round(base_row.mean_n)), ","), int(round(base_row.mean_pos)))]
        for model in ("comp", "cnn", "nt"):
            vals.append("%.3f (%.2f×)" % (base_row[model + "_auprc"], base_row[model + "_lift"]))
        md = macro[(macro.subset == sub) & (macro.dact == dact)].iloc[0]
        vals.append("%+.4f (%+.4f, %+.4f)" % (md.point_delta, md.ci_lo, md.ci_hi))
        for i, (v, w) in enumerate(zip(vals, widths2)):
            cell = cells[i]; cell.width = Cm(w)
            p = cell.paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT if i < 2 else WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_before = Pt(2)
            p.paragraph_format.space_after = Pt(2)
            p.paragraph_format.line_spacing = 1.0
            r = p.add_run(v); r.font.size = Pt(7.5)
            if (sub, dact) == ("unseen", "dev_active"):
                r.bold = True
    for cell in t2.rows[-1].cells:
        set_cell_edge(cell, "bottom", 8)

    # ---------------- Table 3 ----------------
    cap = doc.add_paragraph()
    cap.paragraph_format.page_break_before = False
    cap.paragraph_format.space_before = Pt(6)
    cap.paragraph_format.space_after = Pt(4)
    cap.paragraph_format.line_spacing = 1.15
    cap.paragraph_format.keep_with_next = True
    add_rich(cap, "**Table 3 | Stronger composition and motif baselines.** AUPRC for logistic composition, composition-XGBoost, motif-XGBoost, composition+motif-XGBoost and NT-LoRA in the four strata. Δ is the macro NT − composition-XGBoost difference with its 95% CI.", size=9)
    comp = pd.read_csv(os.path.join(LH, "LH_model_compare_macro.tsv"), sep="\t")
    strong = pd.read_csv(os.path.join(LH, "LH_nt_vs_strong_macro.tsv"), sep="\t")
    boot = pd.read_csv(os.path.join(LH, "LH_nt_vs_strong_boot.tsv"), sep="\t")
    cols3 = ["Locus set", "Development activity", "Composition (logistic)",
             "Composition (XGBoost)", "Motif (XGBoost)",
             "Composition+motif (XGBoost)", "NT-LoRA"]
    widths3 = [2.2, 2.2, 2.4, 2.4, 2.3, 2.6, 2.1]
    t3 = doc.add_table(rows=1, cols=len(cols3))
    t3.alignment = WD_TABLE_ALIGNMENT.CENTER
    t3.autofit = False
    set_table_borders(t3)
    hdr = t3.rows[0]
    trPr = hdr._tr.get_or_add_trPr()
    th = OxmlElement("w:tblHeader"); th.set(qn("w:val"), "true"); trPr.append(th)
    for i, text in enumerate(cols3):
        cell = hdr.cells[i]; cell.width = Cm(widths3[i])
        p = cell.paragraphs[0]; p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(3); p.paragraph_format.space_after = Pt(3)
        p.paragraph_format.line_spacing = 1.0
        r = p.add_run(text); r.bold = True; r.font.size = Pt(7.0)
        r.font.color.rgb = RGBColor(0, 0, 0)
        set_cell_edge(cell, "top", 8)
        set_cell_edge(cell, "bottom", 6)
    for sub, dact in order:
        vals = [SUBSET[sub], DACT[dact]]
        for model in ("comp_logit", "comp_xgb", "motif_xgb", "comb_xgb"):
            v = comp[(comp.subset == sub) & (comp.dact == dact) &
                     (comp.model == model)].auprc.iloc[0]
            vals.append("%.3f" % v)
        srow = strong[(strong.subset == sub) & (strong.dact == dact)].iloc[0]
        vals.append("%.3f" % srow.nt)
        cells = t3.add_row().cells
        for i, (v, w) in enumerate(zip(vals, widths3)):
            cell = cells[i]; cell.width = Cm(w)
            p = cell.paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT if i < 2 else WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_before = Pt(2)
            p.paragraph_format.space_after = Pt(2)
            p.paragraph_format.line_spacing = 1.0
            r = p.add_run(v); r.font.size = Pt(7.0)
            if (sub, dact) == ("unseen", "dev_active"):
                r.bold = True
    for cell in t3.rows[-1].cells:
        set_cell_edge(cell, "bottom", 8)
    parts = []
    for sub, dact in order:
        srow = strong[(strong.subset == sub) & (strong.dact == dact)].iloc[0]
        brow = boot[(boot.subset == sub) & (boot.dact == dact)].iloc[0]
        parts.append("%s %s: %+.4f (%+.4f to %+.4f)" %
                     (SUBSET[sub].lower(), DACT[dact].replace("_", "-"),
                      srow.nt_minus_comp_xgb, brow.comp_xgb_ci_lo, brow.comp_xgb_ci_hi))
    row = t3.add_row()
    cell = row.cells[0].merge(row.cells[-1])
    p = cell.paragraphs[0]
    p.paragraph_format.space_before = Pt(4)
    p.paragraph_format.space_after = Pt(4)
    add_rich(p, "Δ (NT − composition-XGBoost): " + "; ".join(parts), size=7.5)
    set_cell_edge(cell, "bottom", 8)
    move_tables_before_discussion(doc)
    doc.save(OUT)
    print("saved", OUT)


if __name__ == "__main__":
    main()
