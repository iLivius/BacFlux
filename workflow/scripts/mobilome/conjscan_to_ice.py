#!/usr/bin/env python3
"""Turn CONJscan machinery hits + Bakta annotation into ICE / IME candidates.

This is work package E, PHASES 0, 1, 2 and 4 of the mobilome spec
(`docs/mobilome_module_SPEC.md` §8) - what the spec calls the "usable v0":

    Phase 0  data contract   -> every input becomes the same tidy shape
                                (contig, start, end, strand, type, label)
    Phase 1  anchors         -> relaxase / coupling protein / T4SS(MPF) from
                                CONJscan, integrase from the Bakta products
    Phase 2  candidate seeds -> anchors on the same contig, clustered by distance
    Phase 4  classification  -> ICE / IME / island / conjugative region

PHASE 3 (att-site and direct-repeat search, i.e. the real BOUNDARIES of the
element) IS DELIBERATELY NOT IMPLEMENTED HERE. It belongs to the long-read
workflows, because on a short-read assembly the element usually does not sit on
one contig in the first place. Every candidate therefore reports
`boundary_method = none` and empty `attL`/`attR`: the interval below is the span
of the machinery we could see, NOT the true extent of the element. Those three
columns are kept in the schema so the long-read version can fill them in without
changing the file format.

WHY THIS SCRIPT EXISTS (the biology)
------------------------------------
The naive reading of a genome is "chromosome = stays put, plasmid = can move".
That is wrong, and this script is the only place in the module that says so.
An INTEGRATIVE AND CONJUGATIVE ELEMENT (ICE) sits in the chromosome like any
other stretch of DNA, but carries its own conjugation machinery, so it can
excise itself, transfer into another cell and integrate there. A resistance gene
inside an ICE is chromosomal AND transferable at the same time. Plasmid mobility
is answered elsewhere (Platon, plus the plasmid concordance step); the CHROMOSOME
case is answered here and nowhere else.

Three pieces of machinery decide what an element can actually do:

  * RELAXASE (CONJscan `T4SS_MOB*`) - nicks the DNA at the origin of transfer
    and pilots the single strand across. No relaxase, no transfer of that DNA.
  * COUPLING PROTEIN, T4CP (`T4SS_t4cp1/t4cp2`, `T4SS_tcpA`) - the motor that
    hands the relaxase-DNA complex to the secretion system. Necessary, but on
    its own it moves nothing.
  * T4SS / MPF (`T4SS_virb4`, `T4SS_F_tra*`, `T4SS_T_virB*`, ...) - the
    mating-pair formation apparatus, the physical bridge between the two cells.

    ** The coupling protein is NOT counted as a T4SS/MPF component here. **
    That distinction is the whole difference between "can move itself" and
    "can be moved by someone else's machinery" (spec §2.5, ladder tiers 5 vs 6),
    and it is exactly what the real validation sample showed: sample 386 has a
    relaxase plus a coupling protein on the chromosome and NO mating-pair
    apparatus, so it is mobilisable, not self-transmissible.

    ** A mating-pair hit only counts as an APPARATUS when CONJscan called a
    typed T4SS system (`T4SS_type*` / `dCONJ_type*`). ** The `MOB` model is
    relaxase-only by definition and lists VirB4 as an ACCESSORY gene, so one
    incidental VirB4 inside a MOB system is not a mating bridge and must not
    raise an IME to an ICE. The hit is still reported in `has_t4ss`, with the
    reason for not counting it written to the audit file.

  * INTEGRASE (from the Bakta product text) - the recombinase that puts the
    element into the chromosome and takes it out again. It is what separates an
    ICE (integrates) from a plain conjugative region sitting on a contig.

DATA FLOW
---------
Inputs:
  --conjscan-tsv    CONJscan/MacSyFinder 2.1.6 `best_solution.tsv` (rule
                    conjscan). One row per machinery gene found, grouped into
                    systems by `sys_id`, with `sys_wholeness` saying how complete
                    each system is. MAY BE ABSENT - most environmental isolates
                    carry no conjugative system at all, and that is a result, not
                    an error. A directory may be given instead of the file, in
                    which case `best_solution.tsv` is looked for inside it,
                    because that is where MacSyFinder's --out-dir puts it.
  --bakta-gff       Bakta GFF3 for the same sample (rule bakta). Used for TWO
                    things: (a) CONJscan reports protein IDs, not coordinates,
                    so the GFF is what turns a hit into contig/start/end/strand;
                    (b) the integrase anchors are found by matching the CDS
                    `product` text, since CONJscan has no integrase model.
  --contig-lengths  optional `contig<TAB>length` table. When it is not given the
                    lengths are taken from the GFF's own `##sequence-region`
                    lines, which Bakta always writes. Needed only for the
                    contig-end honesty flags.

Outputs:
  --out-table   one row per candidate element, in the column names the sibling
                script `colocalise.py` reads (contig / start / end / strand /
                element_type / mge_id / mge_name), so the co-localisation step
                can consume this file directly as its --is-table alongside the
                ISEScan rows.
  --out-audit   one row per decision: every input that was missing, every hit we
                could not place, every cluster we dropped and every confidence
                downgrade, each with an explicit reason. This is the BacFlux
                convention (`contig_taxonomy_decisions.tsv`): a filtering
                decision that is not written down did not happen.

HOW `element_type` MAPS ONTO THE DOWNSTREAM VOCABULARY (read before changing it)
-------------------------------------------------------------------------------
`colocalise.py` recognises exactly two of our four classes and raises the AMR
gene's mobility tier for them:

    mge_class            element_type written    what colocalise.py does
    ------------------   ---------------------   -----------------------------
    ice                  ice                     tier 6, PREDICTED self-transmissible
    ime                  ime                     tier 5, mobilisable with a helper
    cime_or_island       genomic_island          not recognised -> ignored, and the
                                                 reason is written to ITS audit file
    conjugative_region   conjugative_region      likewise ignored

That is deliberate, not an oversight. A decayed island with only an integrase is
PASSIVE, and a conjugative region with no integrase has no established
boundaries, so neither may be allowed to raise a resistance gene's mobility
tier. They are still written to the table (the spec says "report it, do not call
it an ICE") and the true class is always available in the `mge_class` column.
Giving them an element_type that colocalise.py does not know means they are
skipped LOUDLY, with a line in its audit file, rather than silently.

LANGUAGE DISCIPLINE (spec §2.6)
-------------------------------
Everything here is a PREDICTION from sequence. The output always says "predicted
self-transmissible", never "transmissible"; the confirmatory experiment is a
filter or broth mating assay, not software. Nothing in this file may be labelled
"EFSA-compliant" - it is supporting evidence for the intrinsic/acquired
judgement, not the judgement.

Coordinates are 1-BASED AND INCLUSIVE throughout, which is what GFF3 uses, so no
conversion happens anywhere. Distances are the number of bases strictly BETWEEN
two features, so touching features are 0 bp apart.

Licensing note: an independent implementation. No code was taken from EBI's
mobilome-annotation-pipeline or from ICEfinder (parts of which are CC BY-NC-SA
and incompatible with BacFlux's MIT licence, spec §11). Only conventions are
reused - the `contig|type-start:end` identifier format and the discard-with-
reason file - and conventions are facts, not expression.

Standalone CLI, Python standard library only, exercised by test_conjscan_to_ice.py
with small hand-built fixtures plus the real CONJscan output of sample 386.
"""

import argparse
import csv
import os
import re
import sys

# Phase 3 lives in its own module so the algorithm can be read and tested on
# its own; this script supplies it with the candidates and the sequence.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import att_search
from urllib.parse import unquote


# ── The anchor vocabulary ────────────────────────────────────────────────────
# Four kinds of evidence, and every decision below is made from which of them are
# present. Kept as constants so a typo is a NameError instead of a silent miss.

ANCHOR_RELAXASE = "relaxase"        # T4SS_MOB* - nicks the DNA, pilots the strand
ANCHOR_T4CP = "t4cp"                # coupling protein - the motor, not the bridge
ANCHOR_T4SS = "t4ss"                # mating-pair formation apparatus (MPF)
ANCHOR_INTEGRASE = "integrase"      # site-specific recombinase, from Bakta
# AICE machinery: an actinomycete element that moves as DOUBLE-STRANDED DNA
# between Streptomyces hyphae, pushed through the septal pore by a FtsK/SpoIIIE
# translocase. It has no relaxase and no mating bridge, so it is NOT conjugation
# and this class is never counted as conjugation machinery anywhere below.
ANCHOR_AICE = "aice_machinery"

# The order anchor classes are listed in, so the `anchor_classes` column reads
# the same way for every row and can be compared between samples.
ANCHOR_CLASS_ORDER = [ANCHOR_INTEGRASE, ANCHOR_RELAXASE, ANCHOR_T4CP, ANCHOR_T4SS,
                      ANCHOR_AICE]

# Where one hit came from, so that every downstream rule can ask "was this
# CONJscan or ICEscan?" and the audit trail can say so in words. See the
# ICEscan block further down for why the two are not interchangeable.
SOURCE_CONJSCAN = "conjscan"
SOURCE_ICESCAN = "icescan"
# Not a model set: the Bakta product text, which is where most integrases are
# still found. On the 30 curated pilot elements, 12 have an integrase that ONLY
# the product-text regex finds and no ICEscan profile does, so this is a
# first-class source of evidence and not a fallback.
SOURCE_BAKTA_PRODUCT = "bakta_product"


# ── Which CONJscan gene name is which piece of machinery ─────────────────────
# CONJscan profile names all start with "T4SS_", so the prefix says nothing; the
# REST of the name is what identifies the component. Verified against the
# installed CONJScan 2.1.0 profile list (125 profiles): the relaxases are
# T4SS_MOB{B,C,F,H,M,P1,P2,P3,Q,T,V}, the coupling proteins are T4SS_t4cp1,
# T4SS_t4cp2 and T4SS_tcpA, and everything else (T4SS_virb4, T4SS_F_tra*,
# T4SS_T_virB*, T4SS_G_tfc*, FATA_prg*, ...) is a mating-pair component.
RELAXASE_NAME_PATTERN = re.compile(r"mob", re.IGNORECASE)
T4CP_NAME_PATTERN = re.compile(r"t4cp|tcpa", re.IGNORECASE)


# ── ICEscan's extra vocabulary ───────────────────────────────────────────────
# ICEscan is a FORK of CONJScan 2.0.1 by the same Pasteur authors (its
# metadata.yml still says "CONJScan"), one minor version behind the CONJScan
# 2.1.0 we run. It ADDS an IME model, an AICE model and 21 profiles; it REMOVES
# MOB.xml, the decayed dCONJ models and the whole Plasmids set, and its T4SS
# quorum is stricter. So it is a UNION partner, never a replacement - swapping
# to it loses elements. We run both and take from ICEscan only what it is
# genuinely better at: integrase anchors, the Gram-positive/IME relaxase
# families, and the IME/AICE model classes.
#
# The profiles below are the ones CONJScan 2.1.0 does not have, so
# anchor_class_for_gene_name would otherwise fall through to its default and
# call every one of them a mating-pair component. That default is right for
# CONJscan (all 125 of its profiles really are T4SS_*) and catastrophic here: it
# would invent a mating bridge out of an integrase and promote elements to
# "predicted self-transmissible" on nothing at all.

# INTEGRASE PROFILES WE TRUST.
# These four are genuine element integrases: the tyrosine-recombinase profile
# itself, the large serine recombinase family, and two families ICEscan curated
# for elements CONJScan never modelled.
ICESCAN_TRUSTED_INTEGRASE_PROFILES = {
    "Phage_integrase",   # the tyrosine recombinase of most ICEs and prophages
    "Recombinase",       # large serine recombinase - the Tn916/AICE integrase type
    "UPF0236",           # ICEscan-curated integrase family
    "PB001819",          # ICEscan-curated integrase family
}

# INTEGRASE PROFILES WE REFUSE, each with the protein it really is and what
# trusting it cost us. All four are listed by ICEscan's own IME.xml as
# exchangeables of Phage_integrase, so THIS IS A DELIBERATE DIVERGENCE FROM
# UPSTREAM, not an oversight - a reader comparing us to ICEscan must be able to
# see that we disagreed on purpose, and why.
ICESCAN_EXCLUDED_INTEGRASE_PROFILES = {
    # IntI1, the class-1 INTEGRON integrase. An integron is a gene-capture
    # system, not an integrative element: it has no att site of its own to
    # delimit and it does not excise. Measured on CP042858.1 - an att-bounded
    # 32,103 bp tier-6 ICE became a 103,303 bp UNBOUNDED one at unchanged
    # 'high' confidence, anchored on IntI1. All four of this profile's hits on
    # the 28-genome benchmark are Bakta-annotated "class 1 integron integrase
    # IntI1" or "integron integrase".
    "TIGR02249": "IntI1, the class-1 integron integrase",
    # XerC and XerD, the pair of chromosomal housekeeping recombinases that
    # resolve chromosome dimers at the dif site. Every bacterium has both, they
    # sit wherever the dif site is, and they have nothing to do with mobile
    # elements. Measured: a XerC attached 43,639 bp from a 945 bp relaxase
    # cluster produced an element six times its annotated size. They fire 11 and
    # 9 times raw on this benchmark and are used by nothing, so the damage is
    # real but would be invisible in the final table.
    "TIGR02224": "XerC, a chromosomal dif-site recombinase",
    "TIGR02225": "XerD, a chromosomal dif-site recombinase",
    # The rve integrase catalytic domain - which is the DDE domain shared by IS
    # TRANSPOSASES. Of its 71 hits on this benchmark, 66 are Bakta-annotated
    # transposases (IS6 IS15DIV x14, IS3 family x14, IS6 IS1216E x8, IS30 x6,
    # ...). Note that find_integrase_anchors ALREADY throws these out on the
    # product-text side (TRANSPOSASE_PRODUCT_PATTERN); admitting rve as an HMM
    # integrase would quietly reintroduce exactly the false positives that guard
    # exists to remove, and an IS sitting next to a relaxase would be promoted
    # to an ICE.
    "rve": "the rve/DDE catalytic domain shared by IS transposases",
}

# AICE MACHINERY. An AICE (actinomycete integrative and conjugative element)
# does not conjugate: it encodes a FtsK/SpoIIIE translocase that pushes
# double-stranded DNA through the septal cross-wall into the next hyphal
# compartment, plus its own replication initiator. Neither is conjugation
# machinery, so these are given their own anchor class and are never counted
# as a relaxase or a mating bridge.
#
# WHY THEY MAY NOT SEED AN ELEMENT ON THEIR OWN: FtsK/SpoIIIE is a core
# chromosome-partitioning protein present in essentially every bacterium. A
# cluster is only allowed to become an AICE when ICEscan actually assembled its
# AICE model there (see keep_or_drop_cluster), which requires the translocase,
# a Rep protein and an integrase together.
ICESCAN_AICE_PROFILES = {
    "FtsK_SpoIIIE",   # the translocase that does the moving
    "RepSAv2",        # replication initiator of the pSAM2/SLP1 family
    "Prim-Pol",       # primase-polymerase, exchangeable for RepSAv2
    "DUF3631",        # ditto, exchangeable for RepSAv2
}

# ICEscan's relaxase families, which are what let it see Gram-positive and
# IME relaxases that CONJScan's T4SS_MOB* profiles miss. Every one of them is
# spelled "Relaxase_*", so the prefix is the rule; the two that also contain
# "MOB" (Relaxase_firmi_MOBL, Relaxase_profile_MOBT) and ICEscan's T4SS_MOBL
# are already caught by RELAXASE_NAME_PATTERN and are not special-cased twice.
ICESCAN_RELAXASE_PREFIX = "Relaxase_"


def anchor_class_for_gene_name(gene_name):
    """Say which machinery class one MacSyFinder `gene_name` belongs to.

    Input:  the `gene_name` cell of best_solution.tsv, e.g. 'T4SS_MOBP1' from
            CONJscan or 'Phage_integrase' from ICEscan. NOTE this is the profile
            that ACTUALLY matched, not the model's reference gene - ICEscan
            reports IntI1 hits as gene_name=TIGR02249 under
            hit_gene_ref=Phage_integrase, so reading the reference gene instead
            would let every excluded profile back in through the front door.
    Output: one of the ANCHOR_* constants above, or None for a profile we have
            decided not to trust as evidence of anything (the caller drops it
            and writes an audit row naming what it really is).

    The order of the tests matters.

    The integrase decisions come FIRST, because both lists below name specific
    profiles and must not be reached by any of the pattern tests underneath.
    Then relaxases, because the relaxase is the component that decides whether
    DNA can be transferred at all; then coupling proteins, so that they are NOT
    swept into the T4SS bucket (see the module docstring - that mistake would
    turn every mobilisable element into a self-transmissible one); then the AICE
    translocase, which is not conjugation machinery at all.

    Anything left over is treated as a mating-pair component. That remains the
    safe default because, after the named sets above have been taken out, what
    is left in either tool's vocabulary really is the structural genes of the
    apparatus (T4SS_T_virB*, T4SS_F_tra*, T4SS_G_tfc*, FATA_*, FA_orf*).
    """
    if gene_name in ICESCAN_EXCLUDED_INTEGRASE_PROFILES:
        return None
    if gene_name in ICESCAN_TRUSTED_INTEGRASE_PROFILES:
        return ANCHOR_INTEGRASE
    if gene_name.startswith(ICESCAN_RELAXASE_PREFIX):
        return ANCHOR_RELAXASE
    if RELAXASE_NAME_PATTERN.search(gene_name):
        return ANCHOR_RELAXASE
    if T4CP_NAME_PATTERN.search(gene_name):
        return ANCHOR_T4CP
    if gene_name in ICESCAN_AICE_PROFILES:
        return ANCHOR_AICE
    return ANCHOR_T4SS


def machinery_label(gene_name):
    """Short human label for a machinery hit: 'T4SS_MOBP1' -> 'MOBP1'.

    Only used for the `relaxase_type` / `mpf_type` columns and the audit text,
    so that a report can say "MOBF relaxase" instead of quoting a profile name.
    """
    return gene_name.split("_")[-1] if gene_name else "NA"


# The one alternative name a VirB4 hit can be reported under. Every CONJscan
# model that contains VirB4 lists `T4SS_I_traU` as an EXCHANGEABLE profile for it
# (see MOB.xml, T4SS_typeI.xml and the rest of the definitions directory), so a
# VirB4 can legitimately reach us under either name.
#
# Careful: `T4SS_F_traU` is a DIFFERENT profile - an F-type mating-pair gene in
# its own right, not a VirB4 stand-in - which is why the name below is matched
# exactly rather than as a substring of "traU".
VIRB4_EXCHANGEABLE_GENE_NAME = "T4SS_I_traU"


def hit_is_virb4(hit):
    """Say whether one CONJscan hit is VirB4, under whichever name it was reported.

    Input:  one hit dict from read_conjscan_hits.
    Output: True/False. Used only by the truncation check in
            build_conjscan_anchors, because VirB4 - the ATPase that powers the
            mating bridge - is, with the relaxase, one of the two components that
            has to WORK for conjugation to happen.

    Why this needs two columns rather than one string test: MacSyFinder writes
    the profile that ACTUALLY matched in `gene_name` and keeps the model's own
    gene in `hit_gene_ref`. The real output of sample 006 shows exactly that -
    `gene_name = T4SS_MOBM` for a hit found under the model gene `T4SS_MOBB`. So
    a truncated VirB4 that matched the exchangeable traU profile is written as
    `T4SS_I_traU`, and a check that only looked at `gene_name` for the string
    'virb4' would miss it and leave the element reported as intact.
    """
    gene_name = hit["gene_name"]
    # .get, because `hit_gene_ref` is not one of the columns we insist on: a
    # MacSyFinder version that does not write it must not crash the sample.
    gene_ref = hit.get("hit_gene_ref", "")

    if "virb4" in gene_name.lower() or "virb4" in gene_ref.lower():
        return True

    # Fallback when there is no hit_gene_ref to consult: the exchangeable profile
    # name itself, matched exactly (see the note above about T4SS_F_traU).
    return gene_name.strip() == VIRB4_EXCHANGEABLE_GENE_NAME


def mpf_type_from_model(model_fqn):
    """Pull the MPF type letter out of a CONJscan model name, if it has one.

    Input:  `model_fqn`, e.g. 'CONJScan/Chromosome/T4SS_typeF' or
            'CONJScan/Chromosome/MOB' or 'CONJScan/Chromosome/dCONJ_typeFA'.
    Output: ('F', False) / ('', False) / ('FA', True) - the type letter(s) and
            whether this is a DECAYED system model.

    Two facts from the installed CONJScan 2.1.0 definitions:
      * the `MOB` model is relaxase-centred and has no MPF type at all, so a hit
        under it tells us there is machinery but not which family;
      * the `dCONJ_type*` models describe DECAYED (degenerate) systems. A hit
        under one of those is the tool telling us the machinery is broken, which
        is exactly the overcalling trap the spec warns about (§8 Phase 4), so we
        pass that straight through to the degradation check.
    """
    model_name = model_fqn.split("/")[-1] if model_fqn else ""
    is_decayed = model_name.lower().startswith("dconj")
    if "_type" in model_name:
        return model_name.split("_type", 1)[1], is_decayed
    return "", is_decayed


# The two ICEscan models CONJScan 2.1.0 has no equivalent of. These are the only
# thing we take from ICEscan's system-level output; everything else about a
# system (its span, its wholeness, its mating-pair type) is read from CONJscan
# alone - see build_conjscan_anchors for why.
ICESCAN_MODEL_IME = "IME"
ICESCAN_MODEL_AICE = "AICE"


def icescan_model_class(model_fqn):
    """Return 'IME', 'AICE' or '' for one ICEscan model name.

    Input:  `model_fqn` as MacSyFinder wrote it, e.g. 'ICEscan/Chromosome/IME'
            or 'ICEscan/Chromosome/T4SS_typeF'.
    Output: the bare model name when it is one of the two classes ICEscan adds,
            otherwise '' - its T4SS_type* models duplicate CONJscan's and are
            deliberately not used to classify anything.

    Used by build_candidates to let ICEscan name a class our own anchor-presence
    rules cannot reach (AICE), and to record agreement or disagreement on IME.
    """
    model_name = model_fqn.split("/")[-1] if model_fqn else ""
    if model_name in (ICESCAN_MODEL_IME, ICESCAN_MODEL_AICE):
        return model_name
    return ""


# ── Which Bakta product text counts as an integrase ──────────────────────────
# CONJscan has no integrase model, so the other anchor class comes from the
# annotation text. The pattern is the one the spec fixes at §8 Phase 1.
#
# NOTE what is NOT in it: a bare "recombinase". Bakta annotates RecA as
# "recombinase RecA" and it appears in every genome; RecA is the homologous-
# recombination protein and has nothing to do with site-specific integration, so
# only the "tyrosine recombinase" / "serine recombinase" phrasings match.
#
# WHY THE LIST IS LONGER THAN THE SPEC'S: the spec (§8 Phase 1) gives the pattern
# in terms of protein FAMILY names, but Bakta writes UniRef product text, which
# says the same thing several other ways. Checked against every distinct
# integrase-like product in the KPNIH1 positive control; the two that the
# spec-literal pattern missed are marked below, and missing them is not cosmetic
# — "DNA integration/recombination/inversion protein" is how Bakta annotates
# KPNIH1_04511, the integrase 5.7 kb from KPNIH1's conjugative region. Without it
# that element scored as a bare "conjugative_region" (report, do not call an ICE)
# instead of the ICE it is, so the positive control silently under-called its own
# headline result.
INTEGRASE_PRODUCT_PATTERN = re.compile(
    r"tyr(osine)? recombinase"          # "tyrosine recombinase XerC"; MISSED "Tyr recombinase domain-containing protein"
    r"|phage[_ ]integrase"
    r"|xerc|xerd"
    r"|serine recombinase"
    r"|site.specific recombinase"       # hyphen or space
    r"|dna integration"                 # MISSED "DNA integration/recombination/inversion protein"
    r"|integrase"                       # catches "Integrase", "Integrase family protein", "integron integrase IntI1"
    # ── Domain-name spellings, added 2026-07-28 ──────────────────────────────
    # Bakta often names a tyrosine integrase by its DOMAIN rather than by the
    # word "integrase", and when it does, none of the patterns above fire.
    #
    # This is not hypothetical: it cost us the actual ICEKp on BOTH clinical
    # K. pneumoniae isolates. The integrase of TUM24772's ICEKp is annotated
    # "DUF4102 domain-containing protein" at 1,929,249-1,930,511 - 162 bp from a
    # tRNA-Asn, immediately downstream of the yersiniabactin genes (salicylate
    # synthase Irp9 at 1,927,751-1,929,055), which is exactly where ICEKp's P4-type
    # integrase is documented to sit. With no integrase anchor the element was
    # classified "conjugative_region" ("report it, do not call it an ICE") instead
    # of the ICE it is - so the genuine, publishable finding was demoted, and the
    # ICE label landed on a DIFFERENT element 1.1 Mb away.
    #
    # DUF4102 = Pfam PF13356 = the arm-type DNA-binding domain of phage/P4-family
    # tyrosine recombinases. It is specific to that family, which is why it is safe
    # to match on while a blanket "domain-containing protein" would not be: the
    # same genomes also carry BON, DUF554, HTH-luxR, VTT and Rho
    # "domain-containing" products, none of them integrases.
    r"|duf4102"
    r"|arm[- ]dna[- ]binding"
    r"|arm[- ]type dna[- ]binding",
    re.IGNORECASE,
)

# Known false positive of the pattern above: IS3-family and IS630-family
# TRANSPOSASES are routinely annotated with an "integrase catalytic domain"
# (the rve domain is shared). Counting one as an integrase would let a plain
# insertion sequence sitting next to a relaxase be promoted to an ICE. When the
# product says "transposase" as well, the transposase reading wins and the CDS
# is not used as an anchor; the count is written to the audit file so the
# exclusion is never invisible.
TRANSPOSASE_PRODUCT_PATTERN = re.compile(r"transposase", re.IGNORECASE)



def read_replicon_calls(path):
    """contig -> 'chromosome' | 'plasmid' | 'unknown', from the replicon table.

    WHERE IT COMES FROM: rule mobilome_replicons, which combines Platon's
    chromosome/plasmid split and RDS score with geNomad's independent call via
    the D9 concordance. That is far better evidence than any single marker gene:
    plenty of real plasmids carry no recognisable repA, which is exactly why
    dnaapler reported "No_MMseqs2_hits" for the two Col plasmids on these
    isolates.

    A missing file yields an empty dict and every contig is treated as unknown,
    which leaves behaviour unchanged - the module degrades rather than failing.
    """
    calls = {}
    if not path or not os.path.isfile(path):
        return calls
    with open(path, encoding="utf-8", errors="replace") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if len(lines) < 2:
        return calls
    header = lines[0].split("\t")
    try:
        contig_index = header.index("contig")
        call_index = header.index("replicon")
    except ValueError:
        return calls
    for line in lines[1:]:
        fields = line.split("\t")
        if max(contig_index, call_index) >= len(fields):
            continue
        calls[fields[contig_index].strip()] = fields[call_index].strip().lower()
    return calls


# ── Classification (spec §8 Phase 4) ─────────────────────────────────────────
# The four classes, the element_type each is written as (see the module
# docstring for why the last two are written as something colocalise.py does not
# consume), and the mobility wording. The wording is part of the deliverable:
# "predicted" is never dropped.

MGE_CLASS_ICE = "ice"
MGE_CLASS_IME = "ime"
MGE_CLASS_ISLAND = "cime_or_island"
MGE_CLASS_CONJ_REGION = "conjugative_region"
# The third class ICEscan brings, and the mobility ladder has no rung for it -
# see AICE_MOBILITY below.
MGE_CLASS_AICE = "aice"

ELEMENT_TYPE_FOR_CLASS = {
    MGE_CLASS_ICE: "ice",                        # colocalise.py -> tier 6
    MGE_CLASS_IME: "ime",                        # colocalise.py -> tier 5
    MGE_CLASS_ISLAND: "genomic_island",          # deliberately not consumed
    MGE_CLASS_CONJ_REGION: "conjugative_region",  # deliberately not consumed
    # Recognised by colocalise.py as CONTEXT ONLY: it sets the neighbourhood and
    # caps confidence, and can raise no tier. That is the enforcement point for
    # "an AICE never claims a conjugation tier".
    MGE_CLASS_AICE: "aice",
}

# The mobility sentence for an AICE, and the reason it has no tier.
#
# THE BIOLOGY, because this is the one class the spec's ladder does not cover.
# The ladder's top two rungs are both conjugation: tier 5 is "mobilisable by a
# helper's conjugation machinery", tier 6 is "self-transmissible by its own".
# An AICE has neither a relaxase nor a mating-pair apparatus - it is not a
# conjugative element despite the C in its name. It replicates as a circle and
# is translocated as DOUBLE-STRANDED DNA through the septal pore between
# compartments of a Streptomyces mycelium by its own FtsK/SpoIIIE ATPase. So
# putting it at tier 5 or tier 6 would be a FALSE claim in either direction:
# it is not waiting for a helper, and it cannot mate with another cell.
# It therefore gets no tier at all, and says why.
AICE_MOBILITY = (
    "predicted transferable within the mycelium by FtsK/SpoIIIE translocation "
    "(actinomycete AICE); not on the conjugation mobility ladder"
)
AICE_TIER_REASON = (
    "no tier: the mobility ladder's tiers 5 and 6 are both conjugation, and an "
    "AICE conjugates by neither route - it moves as double-stranded DNA between "
    "hyphal compartments. Reporting it at either tier would be a false claim."
)
# An AICE is not a broken ICE, and `missing_components` must not read as though
# it were: the relaxase and mating bridge are absent BY DEFINITION, not missing.
AICE_MISSING_COMPONENTS = (
    "none expected (an AICE has no relaxase and no mating-pair apparatus)"
)

# What every non-AICE row says in `mobility_tier_reason`. The tier itself is
# assigned by colocalise.py, per AMR gene, from element_type plus the gene's own
# context - this table describes elements, not genes, so it never carries one.
TIER_ASSIGNED_DOWNSTREAM_REASON = (
    "assigned per AMR gene by colocalise.py from element_type and context; this "
    "table describes elements, not genes"
)

# Appended to the mobility wording when the machinery is incomplete. Decayed
# elements are common in real genomes and are where naive tools overcall, so the
# downgrade is visible in the sentence itself, not only in a flag column.
DEGRADED_MOBILITY_SUFFIX = " - machinery incomplete"

# A system whose completeness (`sys_wholeness`, the fraction of the model's genes
# that were actually found) is below this is treated as degraded. Convention,
# not biology: 0.7 is the same threshold the spec uses for a truncated profile
# alignment at §8 Phase 4. Changing it changes how loudly we hedge, never which
# genes were found.
WHOLENESS_INTACT_MIN = 0.7

# The same 0.7, applied to a single hit instead of a whole system: a relaxase or
# VirB4 whose HMM alignment covers less than this much of the profile is a
# fragment, and a fragment does not nick DNA or build a pilus (spec §8 Phase 4,
# "truncation check").
PROFILE_COVERAGE_INTACT_MIN = 0.7

# How many DISTINCT anchor classes a candidate needs before it may be called
# high confidence. Three of the four (say integrase + relaxase + T4SS) means the
# call rests on several independent pieces of evidence rather than one hit.
MIN_ANCHOR_CLASSES_FOR_HIGH = 3


def classify_cluster(has_integrase, has_relaxase, has_mpf_apparatus):
    """Turn "which anchors are present" into a class and a mobility statement.

    This is the whole of spec §8 Phase 4, written as a pure function so it can be
    read - and tested - on its own, without building a cluster first.

    Input:  three booleans, from the anchors of ONE candidate cluster. Note that
            the coupling protein is deliberately not one of them: it is recorded
            in the table but it does not change the class, because a coupling
            protein without a relaxase transfers nothing and a coupling protein
            without an MPF has nothing to hand the DNA to.

            `has_mpf_apparatus` is deliberately STRICTER than "a mating-pair
            component was hit somewhere": it means CONJscan called a typed T4SS
            system. See the block that computes it in build_candidates - an
            accessory VirB4 inside a relaxase-only MOB system does not count.
    Output: (mge_class, mobility_wording).

    The spec's table, verbatim:

        integrase + relaxase + T4SS  -> ICE, self-transmissible. The element can
            excise, build a bridge and integrate in the recipient: chromosomal
            location does NOT mean "cannot move".
        integrase + relaxase, no T4SS -> IME, mobilisable. It can be picked up by
            a co-resident conjugative element's machinery, but cannot transfer on
            its own. This is the commonest real finding.
        integrase, no relaxase       -> CIME / genomic island. It integrates and
            sits there; passive.
        no integrase, relaxase + T4SS -> a conjugative REGION. Machinery is
            present but nothing says it is a discrete integrating element, and
            without an integrase we have no reason to believe there are element
            boundaries at all. The spec is explicit: report it, do not call it
            an ICE.

    Two combinations the spec's table does not list, decided here and documented
    so the fallback is never a surprise:
      * integrase + T4SS but NO relaxase - grouped with the island, because
        without a relaxase the DNA cannot be mobilised whatever else is present.
      * relaxase or T4SS alone, no integrase - grouped with the conjugative
        region, with wording that says which half is missing.
    """
    if has_integrase and has_relaxase and has_mpf_apparatus:
        return MGE_CLASS_ICE, "predicted self-transmissible"

    if has_integrase and has_relaxase:
        return MGE_CLASS_IME, "mobilisable (needs a helper)"

    if has_integrase:
        # Integrase but no relaxase: it can integrate, it cannot be mobilised.
        return MGE_CLASS_ISLAND, "passive"

    if has_relaxase and has_mpf_apparatus:
        return MGE_CLASS_CONJ_REGION, "conjugative region, boundaries not established"

    if has_relaxase:
        # Relaxase (often with a coupling protein) but no mating-pair apparatus
        # and no integrase: the DNA can be mobilised by a helper, but we cannot
        # say what the mobilised unit is.
        return MGE_CLASS_CONJ_REGION, "mobilisable region, boundaries not established"

    # Only mating-pair components: a bridge with nothing to send through it.
    return (MGE_CLASS_CONJ_REGION,
            "conjugation machinery without a relaxase, boundaries not established")


# ── What we write ────────────────────────────────────────────────────────────

# One row per candidate element. The first block of names is fixed by
# colocalise.py's parser (contig / start / end / strand / element_type / mge_id /
# mge_name); everything after it is ours. Deliberately NOT used as column names:
# `type`, `complete`, `completeness`, `family`, `cluster` - colocalise.py reads
# those as ISEScan fields and would misinterpret them.
OUTPUT_COLUMNS = [
    "sample",
    "mge_id",                     # contig|class-start:end  (spec §9 ID format)
    "mge_name",                   # curated name - always NA until a naming DB is wired in
    "element_type",               # what colocalise.py reads: ice | ime | aice | genomic_island | conjugative_region
    "mge_class",                  # what we actually called it (ice|ime|aice|cime_or_island|conjugative_region)
    "mobility",                   # the sentence, always "predicted ..." for ICEs
    # This table describes ELEMENTS, not genes, so it never assigns a mobility
    # tier itself - colocalise.py does that per AMR gene. The column exists so
    # that the one class which can NEVER have a tier says so in words instead of
    # leaving a reader to wonder: an AICE is not on the conjugation ladder at
    # all. Never blank; see AICE_TIER_REASON.
    "mobility_tier",
    "mobility_tier_reason",
    "contig",
    "start",                      # 1-based inclusive, first base of the first anchor
    "end",                        # 1-based inclusive, last base of the last anchor
    "length_bp",                  # end - start + 1: the MACHINERY span, not the element
    "strand",                     # + / - when every anchor agrees, else '.' (mixed)
    "n_anchors",
    "n_anchor_classes",           # 1-4; >=3 is one of the conditions for high confidence
    "anchor_classes",             # comma list, in ANCHOR_CLASS_ORDER
    "has_integrase",
    "has_relaxase",
    "has_t4cp",                   # coupling protein: recorded, but never counted as T4SS
    "has_t4ss",
    "relaxase_type",              # MOBF, MOBP1, ... comma-joined if several
    "mpf_type",                   # F, T, FATA, ... from the CONJscan model name
    "mpf_typed_system",           # TRUE when a full T4SS_type*/dCONJ_type* system was called
    "integrase_products",         # the Bakta product text that matched, ' | '-joined
    "anchor_ids",                 # locus tags with their class, so a reader can look them up
    "conjscan_systems",           # sys_id list - the join key back to best_solution.tsv
    "conjscan_models",            # model_fqn list
    # Which model set each piece of evidence came from, and what ICEscan called
    # here. Present so that a reader comparing this table with ICEscan's own
    # output can see immediately where we agreed and where we did not.
    "evidence_sources",           # conjscan / icescan / bakta_product, comma list
    "icescan_model_class",        # IME | AICE | NA - what ICEscan's models said
    "sys_wholeness_min",          # lowest completeness among the contributing systems
    "machinery_intact",           # TRUE/FALSE - FALSE downgrades the mobility wording
    "degraded_reason",            # why machinery_intact is FALSE, or NA
    # How much this row can be trusted, and why. `system` means MacSyFinder
    # assembled a complete model; `profile_hits_only` means it did NOT, and the
    # row was built from individual profile hits because the machinery was too
    # convincing to discard (see read_conjscan_profile_hits). The latter is
    # always capped at low confidence.
    "evidence_level",             # system | profile_hits_only
    # Which of relaxase / coupling protein / mating-pair apparatus is ABSENT.
    # Always filled in, for both evidence levels, so that a class like
    # "passive island" can be read together with the reason for it.
    "missing_components",
    # Phase 3, the att-site search (att_search.py). When a flanking attL/attR
    # pair is found, `start`/`end` above are REPLACED by the element those
    # boundaries define, because that - not the machinery span - is what actually
    # travels when the element moves, and it is what colocalise.py intersects
    # against the AMR genes to decide the cargo. The machinery span is never lost:
    # it is kept verbatim in machinery_start/machinery_end below.
    "boundary_method",            # tRNA | denovo | none  (how start/end were derived)
    "attL",                       # coordinates of the left repeat, or NA
    "attR",                       # coordinates of the right repeat, or NA
    "att_sequence",               # the repeat itself, so a reader can BLAST it
    "att_length_bp",
    "att_mismatches",             # 0 or 1 between the two copies; >1 is not accepted
    "att_trna",                   # the tRNA the element integrated into, when known
    "machinery_start",            # the CONJscan machinery span, always preserved
    "machinery_end",
    "contig_length",
    "dist_to_contig_start",
    "dist_to_contig_end",
    "dist_to_nearest_contig_end",
    "at_contig_boundary",         # TRUE = the element probably runs off the contig
    "spans_contigs",              # TRUE = this system's hits are on more than one contig
    "confidence",                 # high | medium | low
]

# One row per decision. `action` says what happened to the thing, `reason` is a
# short token to filter on, `detail` is the numbers in plain language.
AUDIT_COLUMNS = [
    "sample",
    "contig",
    "start",
    "end",
    "action",     # input_missing | row_skipped | dropped | kept_flagged
    "reason",     # short machine-readable token
    "detail",     # human-readable explanation, with the numbers behind it
]

# The two thresholds the spec leaves configurable, plus the clustering window.
# All three are CONVENTION, not biology - a real ICE can be 200 kb or 15 kb -
# so the measured span is always reported next to the decision.
DEFAULT_WINDOW_BP = 15000
DEFAULT_MIN_ELEMENT_BP = 8000
# The same floor, for clusters with the IME architecture (integrase + relaxase,
# no mating-pair apparatus). Lower because IME machinery is two genes where an
# ICE's is a twenty-gene operon - see keep_or_drop_cluster for the measurement.
# EMPIRICAL: fitted to a handful of Phase 7 observations, not a published bound.
DEFAULT_MIN_IME_ELEMENT_BP = 2000
DEFAULT_MAX_ELEMENT_BP = 500000

# How close to a contig end counts as "probably truncated". Bigger than the
# equivalent number for an insertion sequence (1 kb vs 100 bp) because an ICE is
# tens of kilobases: if the machinery stops 1 kb from the end of the contig, the
# rest of the element is almost certainly in the missing sequence.
DEFAULT_BOUNDARY_BP = 1000


# ── Small shared helpers ─────────────────────────────────────────────────────

def _warn(message):
    """One warning line on stderr, tagged so it can be grepped out of a Snakemake
    log in which many rules are writing at once."""
    sys.stderr.write("[conjscan_to_ice] WARNING: " + message + "\n")


def _text_or_na(value):
    """Return a stripped value, or the literal 'NA' when it is absent or empty.

    Every BacFlux table writes NA rather than an empty cell, so a missing value
    is visible instead of looking like a formatting slip."""
    if value is None:
        return "NA"
    text = str(value).strip()
    return text if text else "NA"


def _tsv_bool(value):
    """Render yes/no/unknown as TRUE / FALSE / NA.

    TRUE and FALSE (not Python's True/False) because these tables get read in R
    at least as often as in Python, and R turns TRUE/FALSE straight into a
    logical column. None means "we could not tell" and stays NA."""
    if value is None:
        return "NA"
    return "TRUE" if value else "FALSE"


def _reads_true(value):
    """Read a TRUE/FALSE/NA cell back as a boolean.

    The inverse of _tsv_bool, for the passes that re-read rows they already
    wrote. Only the literal TRUE counts: NA means "could not tell", and treating
    that as true would let an unknown quietly become a positive claim.
    """
    return str(value).strip().upper() == "TRUE"


def _reads_tristate(value):
    """Read a TRUE/FALSE/NA cell back as True / False / None.

    The three-way counterpart of _reads_true, and the distinction is not
    pedantic. assess_confidence treats at_contig_boundary=None ("no contig length
    is known, so we cannot tell whether the element runs off the end") as a
    reason to cap confidence, while False ("checked, and it does not") is not.
    Collapsing NA to False therefore turns "we did not know" into "we checked and
    it was fine" - an unknown silently becoming a positive claim, which is the
    one direction this module must never drift in.
    """
    text = str(value).strip().upper()
    if text == "TRUE":
        return True
    if text == "FALSE":
        return False
    return None


def _float_or_none(value):
    """Parse a numeric field, returning None instead of raising.

    Used for sys_wholeness and profile coverage, where "not a number" means the
    completeness check simply cannot be made - which is itself reported - rather
    than a reason to crash the sample."""
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def write_tsv(path, columns, rows):
    """Write rows (list of dicts) as a tab-separated table with a header.

    The header is ALWAYS written, even for zero rows. An explicit empty table is
    a result - "this genome has no conjugation machinery", which is the common
    case for environmental isolates - whereas a missing or headerless file looks
    to the next rule like a crash."""
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def audit_row(sample, action, reason, detail, contig="NA", start="NA", end="NA"):
    """Build one audit record. Kept as a function so every caller writes the
    same seven fields in the same order and nothing is forgotten."""
    return {
        "sample": sample,
        "contig": _text_or_na(contig),
        "start": _text_or_na(start),
        "end": _text_or_na(end),
        "action": action,
        "reason": reason,
        "detail": detail,
    }


# ── Phase 0, input 1: the Bakta annotation ───────────────────────────────────

def parse_gff_attributes(attribute_field):
    """Split a GFF3 column-9 attribute string into a plain dict.

    GFF3 percent-encodes any character that would break the format, so Bakta
    writes a product containing a comma as
    `product=Type IV secretory pathway%2C VirD4 component`. We decode it here
    (stdlib `unquote`), because the integrase search matches on the product TEXT
    and an encoded string would not match the way a human reading the file
    expects it to.
    """
    attributes = {}
    for chunk in attribute_field.strip().split(";"):
        if not chunk or "=" not in chunk:
            continue
        key, _, value = chunk.partition("=")
        attributes[unquote(key.strip())] = unquote(value.strip())
    return attributes


def parse_bakta_gff(path):
    """Read a Bakta GFF3 into (features_by_id, cds_features, contig_lengths).

    Input: `05.annotation/bakta/{sample}/{sample}.gff3` (the rule `bakta`).

    What we take from it:
      * every CDS, indexed under BOTH its `ID` and its `locus_tag`. CONJscan
        reports the protein identifier from the .faa it was given, and in a real
        Bakta GFF3 those two attributes carry the same value, e.g.

            contig_2  Pyrodigal  CDS  37327  38967  .  +  0
            ID=386_00040;Name=Relaxase/mobilization nuclease family protein;
            locus_tag=386_00040;product=Relaxase/mobilization nuclease family protein

        so CONJscan's `hit_id` 386_00040 finds contig_2:37327-38967(+). Indexing
        both keys means a future Bakta that changes which one it puts in the
        FASTA header still joins.
      * the CDS list in file order, which the integrase search walks.
      * the contig lengths from the `##sequence-region contig_1 1 4177018`
        header lines, so the contig-end flags work even when no separate
        contig-length table was supplied.

    Parsing stops at `##FASTA`: Bakta appends the whole assembly to the end of
    the GFF3, and those lines are sequence, not features.

    Only CDS rows are indexed. CONJscan hits are proteins, and the integrase
    anchors are proteins too, so tRNA/rRNA/oriC rows cannot contribute an anchor.
    """
    features_by_id = {}
    cds_features = []
    contig_lengths = {}

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break

            if line.startswith("##sequence-region"):
                # '##sequence-region contig_1 1 4177018' -> contig_1 is 4177018 bp
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        contig_lengths[parts[1]] = int(parts[3])
                    except ValueError:
                        pass
                continue

            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                # A feature line without usable coordinates is unusable; skipping
                # it is safer than inventing a position for a machinery gene.
                continue

            attributes = parse_gff_attributes(fields[8])
            feature = {
                "contig": fields[0],
                "start": start,
                "end": end,
                "strand": fields[6] if fields[6] in {"+", "-"} else ".",
                "product": attributes.get("product", ""),
                "feature_id": attributes.get("ID") or attributes.get("locus_tag") or "",
            }
            cds_features.append(feature)

            for key in ("ID", "locus_tag"):
                identifier = attributes.get(key)
                # setdefault: if two CDS somehow share an identifier, the first
                # one wins rather than the last, so the result is stable.
                if identifier:
                    features_by_id.setdefault(identifier, feature)

    return features_by_id, cds_features, contig_lengths


def read_contig_lengths(path):
    """Read an optional `contig<TAB>length` table into {contig: length}.

    Input: whatever the rule passes as --contig-lengths - normally the table the
    ISEScan loader already writes for this assembly, so both scripts measure
    distances against exactly the same contigs. A samtools `.fai` also works,
    because only the first two columns are read.

    A header line is tolerated: any line whose second field is not a number is
    skipped. Returns {} for a path that is missing or empty, and the caller then
    falls back to the GFF's own sequence-region lengths.
    """
    lengths = {}
    if not path or not os.path.exists(path):
        return lengths

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            try:
                lengths[fields[0].strip()] = int(fields[1].strip())
            except ValueError:
                continue          # header row, or a malformed line
    return lengths


# ── Phase 0, input 2: the CONJscan result ────────────────────────────────────

# The columns we actually need out of best_solution.tsv. Matched BY NAME, so
# MacSyFinder adding or reordering a column cannot silently shift our values.
CONJSCAN_REQUIRED_COLUMNS = ["hit_id", "gene_name", "model_fqn", "sys_id"]


def resolve_conjscan_path(path):
    """Accept either the best_solution.tsv itself or MacSyFinder's --out-dir.

    MacSyFinder writes several files into one output directory; the rule may
    hand us the directory. Returns the path to use, or None when there is
    nothing to read - which is a normal result, not an error.
    """
    if not path:
        return None
    if os.path.isdir(path):
        candidate = os.path.join(path, "best_solution.tsv")
        return candidate if os.path.exists(candidate) else None
    return path if os.path.exists(path) else None


def namespaced_sys_id(sys_id, source):
    """Prefix an ICEscan system id so it can never be confused with a CONJscan one.

    WHY THIS IS NOT COSMETIC. MacSyFinder builds a system id as
    {replicon}_{model}_{n}, and the two tools share model names - both define
    T4SS_typeF. Run over the same genome they therefore produce genuinely
    DIFFERENT systems under IDENTICAL ids, e.g. 'CP042858_1_T4SS_typeF_3' from
    each. Merging them would pool their contigs (a false spans_contigs flag),
    pool their wholeness, and let merge_clusters_sharing_a_system join two
    clusters that no single tool ever said belonged together.

    CONJscan ids are left exactly as they are, so nothing about the existing
    single-tool path changes.
    """
    if source != SOURCE_ICESCAN or not sys_id:
        return sys_id
    return f"{SOURCE_ICESCAN}:{sys_id}"


def read_conjscan_hits(path, source=SOURCE_CONJSCAN):
    """Read CONJscan's best_solution.tsv into one dict per machinery hit.

    Input: the file written by `macsyfinder --models CONJScan/Chromosome all`.
    Its real shape (verified on sample 386, MacSyFinder 2.1.6 / CONJScan 2.1.0):
    three '#' banner lines, then a 22-column header, then one row per gene, with
    BLANK LINES separating systems.

    Output: a list of dicts, one per row. Nothing is interpreted here beyond
    parsing - the machinery classification happens in build_conjscan_anchors.

    Returns [] when the file is absent, empty, header-only, or contains only
    comments. All four mean the same biological thing: no conjugation machinery
    was found, which is the common case and must never be an error.

    Raises ValueError only if there ARE data rows but the header is not the one
    we know. That is a tool-version change, and quietly reporting "no ICE" from a
    misparsed file would be much worse than stopping.
    """
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []

    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = [row for row in reader if row and row[0].strip()
                and not row[0].startswith("#")]

    if not rows:
        return []

    header = [name.strip() for name in rows[0]]
    index_of = {name: position for position, name in enumerate(header)}
    missing = [name for name in CONJSCAN_REQUIRED_COLUMNS if name not in index_of]
    if missing:
        raise ValueError(
            f"CONJscan table '{path}' is missing required column(s): "
            f"{', '.join(missing)}. Header was: {header}"
        )

    def field(row, name):
        position = index_of.get(name)
        if position is None or position >= len(row):
            return ""
        return row[position].strip()

    hits = []
    for row in rows[1:]:
        hits.append({
            "hit_id": field(row, "hit_id"),
            "gene_name": field(row, "gene_name"),
            "model_fqn": field(row, "model_fqn"),
            "sys_id": namespaced_sys_id(field(row, "sys_id"), source),
            "sys_wholeness": _float_or_none(field(row, "sys_wholeness")),
            "hit_profile_cov": _float_or_none(field(row, "hit_profile_cov")),
            "hit_status": field(row, "hit_status"),
            # The model's OWN gene, which may differ from gene_name: see
            # hit_is_virb4 below for why we need both.
            "hit_gene_ref": field(row, "hit_gene_ref"),
            # Which model set found this. Carried on every hit so the audit can
            # say where a piece of evidence came from, and so the rules below
            # can treat ICEscan's system-level claims differently.
            "source": source,
        })
    return hits


# ── The fallback: profile hits when MacSyFinder assembled no system ──────────

# MacSyFinder's per-profile hit tables, one file per profile, written next to
# best_solution.tsv. The columns are fixed by the tool; recorded here because the
# parser is written against them.
#
#   # hit_id  replicon_name  position_hit  hit_sequence_length  gene_name
#   i_eval  score  profile_coverage  sequence_coverage  begin  end
HMMER_EXTRACT_SUFFIX = ".res_hmm_extract"


def default_hmmer_dir(conjscan_path):
    """Where MacSyFinder puts hmmer_results/, given best_solution.tsv's path.

    The two always sit side by side in one output directory, so the caller does
    not have to pass both. Returns '' when there is nothing to derive from, and
    read_conjscan_profile_hits then simply finds no hits.
    """
    if not conjscan_path:
        return ""
    base = (conjscan_path if os.path.isdir(conjscan_path)
            else os.path.dirname(conjscan_path))
    return os.path.join(base, "hmmer_results") if base else "hmmer_results"


def read_conjscan_profile_hits(hmmer_dir, source=SOURCE_CONJSCAN):
    """Read the individual profile hits CONJscan found, ignoring system assembly.

    WHY THIS EXISTS. MacSyFinder only reports a SYSTEM when a model's quorum is
    met, and every conjugative model requires a relaxase. On the Phase 7
    benchmark, ICEVflInd1 - a genuine 114 kb SXT/R391-family ICE - hit twelve
    mating-pair profiles (traE, traF, traH, traK, traL, traN, traU, traV, traW),
    a coupling protein AND VirB4, the signature ATPase, but no MOB* profile. No
    system was assembled, so the module reported NOTHING: a real conjugative
    element vanished because one component of the quorum was missing.

    Reporting nothing there is defensible but not useful. Reporting it as though
    a system had been found would be dishonest. So these hits are read as a
    SECOND, WEAKER evidence tier: the element is surfaced, the missing component
    is named in `missing_components`, `evidence_level` says the system was never
    assembled, and confidence is capped at low. The reader can then see both what
    was found and what was not.

    Input:  the hmmer_results/ directory inside CONJscan's output.
    Does:   reads every {profile}.res_hmm_extract. These are ALREADY filtered by
            MacSyFinder's own i-evalue (0.001) and profile-coverage (0.5)
            thresholds, so they are vetted hits rather than raw HMM noise - we
            are not lowering a statistical bar, only skipping the quorum rule.
    Output: hit dicts in exactly the shape read_conjscan_hits returns, so
            everything downstream is unchanged. sys_id is deliberately EMPTY:
            there is no system, and an empty sys_id also stops
            merge_clusters_sharing_a_system from grouping on absent information.

    Returns [] when the directory is missing or holds no hits, which - as
    everywhere else in this module - means "no machinery", never an error.
    """
    hits = []
    if not hmmer_dir or not os.path.isdir(hmmer_dir):
        return hits

    for filename in sorted(os.listdir(hmmer_dir)):
        if not filename.endswith(HMMER_EXTRACT_SUFFIX):
            continue
        path = os.path.join(hmmer_dir, filename)
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as handle:
                lines = [line.rstrip("\n") for line in handle]
        except OSError:
            continue

        for line in lines:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            hits.append({
                "hit_id": fields[0].strip(),
                "gene_name": fields[4].strip(),
                # No model was satisfied, so there is no model_fqn or sys_id to
                # report. Both stay empty and the columns read NA downstream.
                "model_fqn": "",
                "sys_id": "",
                "sys_wholeness": None,
                "hit_profile_cov": _float_or_none(fields[7].strip()),
                "hit_status": "",
                "hit_gene_ref": fields[4].strip(),
                "source": source,
            })
    return hits


# The components a complete conjugative system needs, in the order a reader
# wants them. Used only to say which one was ABSENT when the fallback fires.
CANONICAL_COMPONENTS = [
    (ANCHOR_RELAXASE, "relaxase"),
    (ANCHOR_T4CP, "coupling protein"),
    (ANCHOR_T4SS, "mating-pair apparatus"),
]


def missing_components_for(cluster):
    """Name the canonical machinery components this cluster does NOT have.

    Returned as a comma-joined string for the `missing_components` column, or
    'none' when all three are present. This is the transparency half of the
    fallback: a call built from profile hits alone must say what is missing, so
    that "passive island" can be read as "no relaxase was detected" rather than
    as a positive statement about the element's biology.
    """
    absent = [label for anchor_class, label in CANONICAL_COMPONENTS
              if not cluster_has_class(cluster, anchor_class)]
    return ", ".join(absent) if absent else "none"


# ── Phase 1: anchors ─────────────────────────────────────────────────────────

def make_anchor(feature, anchor_class, label, hit=None):
    """Put one piece of evidence into the module's common shape (spec §8 Phase 0).

    Every anchor - whether it came from CONJscan or from a Bakta product string -
    looks the same from here on: where it is, what class of evidence it is, and
    which CONJscan system (if any) it belongs to. Phase 2 clusters these and
    Phase 4 classifies the clusters.
    """
    return {
        "contig": feature["contig"],
        "start": feature["start"],
        "end": feature["end"],
        "strand": feature["strand"],
        "anchor_class": anchor_class,
        "label": label,
        "feature_id": feature["feature_id"],
        "sys_id": hit["sys_id"] if hit else "",
        "model_fqn": hit["model_fqn"] if hit else "",
        "gene_name": hit["gene_name"] if hit else "",
        # Where this evidence came from: 'conjscan', 'icescan', or - when there
        # is no HMM hit behind it at all - the Bakta annotation text. Read by
        # the integrase tie-break and reported in the audit trail.
        "source": hit.get("source", SOURCE_CONJSCAN) if hit else SOURCE_BAKTA_PRODUCT,
    }


def build_conjscan_anchors(sample, hits, features_by_id):
    """Place every CONJscan hit on the assembly and label what it is.

    Input:  the hit list from read_conjscan_hits, and the CDS index from the
            Bakta GFF3.
    Does:   looks each `hit_id` up in the annotation to get its coordinates, then
            asks anchor_class_for_gene_name what machinery it is. At the same
            time it summarises each SYSTEM (`sys_id`), because completeness and
            the contig-spanning check are properties of the system, not of one
            gene.
    Output: (anchors, integrase_candidates, systems, audit_rows).

    `integrase_candidates` is returned SEPARATELY from the machinery anchors,
    and only ICEscan puts anything in it (CONJScan 2.1.0 has no integrase model
    at all). They are kept apart for two reasons, both of which matter:

      * an integrase belongs at the ELEMENT BOUNDARY, tens of kb from the
        machinery operon, so it needs the wide --integrase-window-bp radius that
        attach_nearby_integrases applies - not the tight machinery window; and
      * ICEscan must not be allowed to set an element's extent. Its IME
        "systems" are not loci - both mandatory genes are declared loners, so
        MacSyFinder exempts them from clustering, and on this benchmark 17 of 71
        of them span more than 1,000 genes (the largest being an entire 10.1 Mb
        chromosome called as one "IME"). Letting such a hit join a machinery
        cluster directly would drag the cluster across the replicon.

    Putting them in the same pool as the Bakta product-text integrases means
    both sources compete under one rule, which is where the FIX-4 tie-break in
    attach_nearby_integrases does its work.

    `systems` is keyed by sys_id and holds:
      wholeness    lowest sys_wholeness seen for it (all rows of a system carry
                   the same value, but taking the minimum is robust)
      decayed      TRUE when the model is one of CONJscan's dCONJ_type* decayed
                   system definitions
      truncated    TRUE when a relaxase or VirB4 hit covers less than 70% of its
                   HMM profile - a fragment of a relaxase nicks nothing
      contigs      the set of contigs its hits landed on. MacSyFinder is run with
                   --db-type ordered_replicon, which treats a whole draft
                   assembly as ONE pseudo-replicon, so it can and does join genes
                   across contig boundaries. A system on more than one contig is
                   therefore an assembly artefact risk and is capped at low
                   confidence later.

    A hit whose protein ID is not in the GFF gets an audit row and is dropped -
    without coordinates it cannot be clustered with anything.
    """
    anchors = []
    integrase_candidates = []
    systems = {}
    audit_rows = []
    # Counted, not listed one by one: `rve` alone fires 71 times on the
    # benchmark, and a row per hit would drown the audit file.
    excluded_profile_counts = {}
    unused_icescan_apparatus = 0

    for hit in hits:
        source = hit.get("source", SOURCE_CONJSCAN)
        system = systems.setdefault(hit["sys_id"], {
            "wholeness": None,
            "decayed": False,
            "truncated": False,
            "contigs": set(),
            "models": set(),
            "mpf_types": set(),
            "n_hits": 0,
            # Which tool assembled this system. build_candidates reads every
            # OTHER field here for CONJscan systems only - see the note there.
            "source": source,
            # 'IME' / 'AICE' / '' - the one thing we do take from an ICEscan
            # system, because CONJScan 2.1.0 models neither class.
            "icescan_class": icescan_model_class(hit["model_fqn"]),
        })
        system["n_hits"] += 1
        system["models"].add(hit["model_fqn"])

        mpf_type, is_decayed = mpf_type_from_model(hit["model_fqn"])
        if mpf_type:
            system["mpf_types"].add(mpf_type)
        if is_decayed:
            system["decayed"] = True

        wholeness = hit["sys_wholeness"]
        if wholeness is not None:
            if system["wholeness"] is None or wholeness < system["wholeness"]:
                system["wholeness"] = wholeness

        anchor_class = anchor_class_for_gene_name(hit["gene_name"])

        # A profile we have decided not to trust (IntI1, XerC, XerD, rve - see
        # ICESCAN_EXCLUDED_INTEGRASE_PROFILES for what each really is and what
        # trusting it cost). It becomes no anchor of any kind: not an integrase,
        # and not a mating-pair component either, which is what the fall-through
        # default would otherwise have made of it.
        if anchor_class is None:
            excluded_profile_counts[hit["gene_name"]] = (
                excluded_profile_counts.get(hit["gene_name"], 0) + 1
            )
            continue

        # ICEscan's mating-pair and coupling-protein hits are READ (they still
        # count towards the system summary above) but are never turned into
        # anchors. Two reasons, and both are deliberate limits on how far we
        # trust the fork:
        #   * its T4SS quorum is stricter than CONJScan 2.1.0's and its model set
        #     drops MOB.xml and the decayed dCONJ models, so its apparatus calls
        #     are a different instrument from the one every threshold here was
        #     tuned against - mixing the two would silently change what counts as
        #     "self-transmissible"; and
        #   * an anchor extends a cluster, and ICEscan system spans are not
        #     element spans (see the docstring above).
        # CONJscan already covers the apparatus better, so nothing is lost.
        if source == SOURCE_ICESCAN and anchor_class in (ANCHOR_T4SS, ANCHOR_T4CP):
            unused_icescan_apparatus += 1
            continue

        # Truncation check (spec §8 Phase 4): only the two components that have
        # to WORK for transfer are tested - the relaxase, and VirB4, the ATPase
        # that powers the mating bridge. A short alignment to a structural
        # accessory gene is much less informative.
        # hit_is_virb4 checks both names VirB4 can arrive under; see the note
        # there about the exchangeable T4SS_I_traU profile.
        is_core_component = (
            anchor_class == ANCHOR_RELAXASE
            or hit_is_virb4(hit)
        )
        coverage = hit["hit_profile_cov"]
        if is_core_component and coverage is not None and coverage < PROFILE_COVERAGE_INTACT_MIN:
            system["truncated"] = True

        feature = features_by_id.get(hit["hit_id"])
        if feature is None:
            audit_rows.append(audit_row(
                sample, "row_skipped", "hit_id_not_in_annotation",
                f"CONJscan reported {hit['gene_name']} for protein "
                f"'{hit['hit_id']}' (system {hit['sys_id']}), but that ID is not a "
                "CDS in the Bakta GFF3, so it has no coordinates and cannot be "
                "placed. CONJscan and Bakta were probably run on different "
                "annotations of this sample.",
            ))
            continue

        system["contigs"].add(feature["contig"])

        # AICE machinery keeps its full profile name as the label ('FtsK_SpoIIIE'
        # rather than machinery_label's 'SpoIIIE'), because the reader needs to
        # recognise the translocase by name.
        label = (hit["gene_name"] if anchor_class == ANCHOR_AICE
                 else machinery_label(hit["gene_name"]))
        anchor = make_anchor(feature, anchor_class, label, hit=hit)

        if anchor_class == ANCHOR_INTEGRASE:
            integrase_candidates.append(anchor)
        else:
            anchors.append(anchor)

    if excluded_profile_counts:
        listed = ", ".join(
            f"{name} x{count} ({ICESCAN_EXCLUDED_INTEGRASE_PROFILES[name]})"
            for name, count in sorted(excluded_profile_counts.items())
        )
        audit_rows.append(audit_row(
            sample, "row_skipped", "untrusted_integrase_profile",
            f"{sum(excluded_profile_counts.values())} hit(s) to profiles that "
            f"ICEscan's IME model lists as integrases but that we do not trust: "
            f"{listed}. We diverge from upstream here ON PURPOSE - none of these "
            "is an element integrase, and each has been measured to inflate or "
            "misplace an element. They were used as no anchor at all.",
        ))

    if unused_icescan_apparatus:
        audit_rows.append(audit_row(
            sample, "row_skipped", "icescan_apparatus_hit_not_used_as_anchor",
            f"{unused_icescan_apparatus} ICEscan mating-pair/coupling-protein "
            "hit(s) were read but not used as anchors. The mating-pair apparatus "
            "is taken from CONJScan 2.1.0 alone: ICEscan is a fork of CONJScan "
            "2.0.1 with a stricter T4SS quorum and no MOB or decayed-system "
            "models, and its system spans are not element spans. ICEscan is used "
            "here for integrase anchors, its Gram-positive relaxase families and "
            "its IME/AICE model classes only.",
        ))

    return anchors, integrase_candidates, systems, audit_rows


def find_integrase_anchors(sample, cds_features):
    """Find the integrase anchors in the Bakta annotation (spec §8 Phase 1).

    Input:  every CDS from the GFF3, in file order.
    Does:   matches the CDS product text against INTEGRASE_PRODUCT_PATTERN, and
            throws out the hits that are really transposases (see the note next
            to TRANSPOSASE_PRODUCT_PATTERN).
    Output: (anchors, audit_rows).

    Why the integrase matters: it is the difference between "there is conjugation
    machinery somewhere on this contig" and "there is an element that can put
    itself into, and take itself out of, a chromosome". Without it we refuse to
    say ICE.

    A caveat worth stating: this is a text match on an annotation, not an HMM
    search. It will find integrases that are pseudogenes, and it will miss ones
    Bakta annotated as "hypothetical protein". Both directions are visible - the
    matched product text is carried into the output table so the reader can
    judge each call.
    """
    anchors = []
    excluded_products = []

    for feature in cds_features:
        product = feature["product"]
        if not product or not INTEGRASE_PRODUCT_PATTERN.search(product):
            continue
        if TRANSPOSASE_PRODUCT_PATTERN.search(product):
            excluded_products.append(product)
            continue
        anchors.append(make_anchor(feature, ANCHOR_INTEGRASE, product))

    audit_rows = []
    if excluded_products:
        # One summary row, not one per gene: a genome can carry dozens of IS
        # transposases and that would drown the audit file.
        audit_rows.append(audit_row(
            sample, "row_skipped", "integrase_match_is_a_transposase",
            f"{len(excluded_products)} CDS matched the integrase product pattern but "
            "also say 'transposase' (IS3/IS630-family transposases share the rve "
            "integrase catalytic domain). They were NOT used as integrase anchors, "
            "because counting one would promote an insertion sequence next to a "
            f"relaxase into an ICE. Examples: "
            f"{' | '.join(sorted(set(excluded_products))[:3])}",
        ))

    return anchors, audit_rows


# ── Phase 2: candidate seeding ───────────────────────────────────────────────


# How far from a machinery cluster an INTEGRASE may sit and still belong to the
# same element. Deliberately much larger than DEFAULT_WINDOW_BP, because the two
# anchor types have different geometry:
#
#   * conjugation machinery is an OPERON - relaxase, coupling protein and the
#     mating-pair genes sit within a few kb of each other, so 15 kb clusters them
#     correctly and a wider window would wrongly merge separate systems;
#   * the integrase sits at the ELEMENT BOUNDARY, next to the tRNA the element
#     integrated into. For a 50-100 kb ICE that is tens of kb from the machinery.
#
# Measured: on BOTH clinical K. pneumoniae isolates the ICEKp integrase (a
# DUF4102/P4-family protein abutting a tRNA-Asn) sits 33,361 bp from the
# machinery cluster. At a shared 15 kb window it was never attached, so the
# element came out as "conjugative_region" - relaxase and T4SS but no integrase,
# i.e. "report it, do not call it an ICE" - and the genuine ICEKp was missed on
# both genomes while the ICE label landed on a different element 1.1 Mb away.
#
# This is a known failure mode, not a local quirk: icefinder-opt's README records
# the same fix upstream - "the tRNA search window was too limited to reliably
# identify ICE boundaries in many cases", widened in their 2025 release.
DEFAULT_INTEGRASE_WINDOW_BP = 50000


# How the sources rank when two integrase candidates are otherwise equal.
# Lower wins. The Bakta product text is preferred over an HMM profile hit at
# equal distance and equal att support, because it names the protein in words a
# reader can check ("class 1 integron integrase IntI1" is a sentence you can
# disagree with; "TIGR02249" is not) and because it is the source that found the
# integrase for 12 of the 30 curated pilot elements.
INTEGRASE_SOURCE_RANK = {
    SOURCE_BAKTA_PRODUCT: 0,
    SOURCE_ICESCAN: 1,
    SOURCE_CONJSCAN: 2,
}


def attach_nearby_integrases(clusters, integrase_anchors, integrase_window_bp,
                             att_supports=None):
    """Give each machinery cluster the integrase that belongs to it, if any.

    Takes in: clusters built from MACHINERY anchors only; every integrase
              candidate from Phase 1 (Bakta product text AND, when ICEscan ran,
              its trusted integrase profiles); and optionally `att_supports`.
    Does:     for each cluster, ranks every integrase on the same contig within
              integrase_window_bp and attaches the single best one.
    Returns:  (clusters, audit_rows) - clusters modified in place.

    Only ONE integrase is attached per cluster. An element has one integrase;
    attaching every integrase within 50 kb would manufacture anchor classes on a
    chromosome that is full of prophage integrases.

    HOW THE BEST ONE IS CHOSEN (this is FIX 4, and it is a correctness fix).
    The old rule was "closest, and first-seen wins a tie". That is fine while
    there is only one source of integrases, but once ICEscan's profile hits join
    the pool the candidates routinely sit INSIDE the machinery cluster, where
    every one of them has a gap of 0. A tie among all candidates meant the choice
    fell to whichever happened to be first in the list - i.e. to the SOURCE
    rather than to any evidence. That is exactly how the CP042858.1 regression
    happened: an att-bounded 32,103 bp ICE became a 103,303 bp unbounded one, at
    unchanged 'high' confidence, because an IntI1 hit won a tie it should never
    have been in.

    So the ranking key is, in order:

      1. does this integrase YIELD AN att PAIR?  An att site is the scar left
         when the element recombined into the chromosome, so an integrase whose
         inclusion produces flanking attL/attR is the one that actually put this
         element here. This is real evidence about which integrase belongs to the
         element, and it is the only key that is - the other two are tie-breaks.
      2. distance to the cluster.
      3. the source, per INTEGRASE_SOURCE_RANK above.
      4. position, purely so the result is deterministic when all else is equal.

    `att_supports(contig, start, end) -> bool` is supplied by the caller when a
    genome FASTA is available; without one it is None, key 1 is constant, and the
    ranking degrades to "closest, then source" - which is the old behaviour made
    deterministic. Note the consequence of key 1 being primary: an att-supporting
    integrase BEATS a closer one that yields nothing. That is intended. Distance
    is a proxy for belonging; an att pair is direct evidence of it.
    """
    audit_rows = []
    attached_ids = set()

    def gap_to(cluster_lo, cluster_hi, anchor):
        """Bases strictly between the cluster and this integrase; 0 if they overlap."""
        if anchor["start"] > cluster_hi:
            return anchor["start"] - cluster_hi
        if anchor["end"] < cluster_lo:
            return cluster_lo - anchor["end"]
        return 0

    for cluster in clusters:
        contig = cluster[0]["contig"]
        lo = min(a["start"] for a in cluster)
        hi = max(a["end"] for a in cluster)

        # Every integrase near enough to be in the running, with its distance.
        in_range = []
        for anchor in integrase_anchors:
            if anchor["contig"] != contig:
                continue
            gap = gap_to(lo, hi, anchor)
            if gap <= integrase_window_bp:
                in_range.append((gap, anchor))
        if not in_range:
            continue

        # Key 1 is the expensive one (it runs the att search), so it is only
        # computed when there is actually a choice to make. With a single
        # candidate the answer cannot change the outcome.
        def att_supported(anchor):
            if att_supports is None or len(in_range) < 2:
                return False
            span_lo = min(lo, anchor["start"])
            span_hi = max(hi, anchor["end"])
            return bool(att_supports(contig, span_lo, span_hi))

        ranked = sorted(
            in_range,
            key=lambda pair: (
                not att_supported(pair[1]),                 # False (0) sorts first
                pair[0],                                    # then closest
                INTEGRASE_SOURCE_RANK.get(pair[1]["source"], 9),
                pair[1]["start"],
            ),
        )
        gap, anchor = ranked[0]

        # Clusters are built from MACHINERY anchors only - integrases, from
        # either source, are held back for exactly this step - so a cluster
        # never already contains one and there is nothing to de-duplicate here.
        cluster.append(anchor)
        attached_ids.add(id(anchor))
        runners_up = len(in_range) - 1
        audit_rows.append(audit_row(
            "", "evidence_recorded", "integrase_attached_beyond_cluster_window",
            f"{contig}:{lo}-{hi}: an integrase at {anchor['start']}-{anchor['end']} "
            f"({anchor.get('label','')}, from {anchor.get('source','')}) sits "
            f"{gap} bp away - beyond the "
            f"machinery clustering window but within the {integrase_window_bp} bp "
            "integrase window. Attached: conjugation machinery is an operon and "
            "clusters tightly, whereas the integrase sits at the element "
            "boundary, tens of kb away on a large ICE."
            + (f" Chosen over {runners_up} other candidate(s) by att support, "
               "then distance, then source." if runners_up else ""),
            contig=contig, start=lo, end=hi))

    # Integrases belonging to NO machinery cluster. Under the old shared-window
    # design these formed integrase-only clusters that were then dropped with a
    # stated reason; clustering machinery alone means they never form one, so the
    # reason has to be recorded here or the decision becomes invisible - and every
    # filtering decision in this module carries a reason.
    orphans = [a for a in integrase_anchors if id(a) not in attached_ids]
    if orphans:
        audit_rows.append(audit_row(
            "", "not_applicable", "integrase_without_conjugation_machinery",
            f"{len(orphans)} integrase(s) were found with no conjugation "
            f"machinery within {integrase_window_bp} bp, so they anchor no "
            "element. An integrase alone is not an ICE - chromosomes carry many, "
            "mostly from prophages. Examples: "
            + "; ".join(f"{a['contig']}:{a['start']}-{a['end']}" for a in orphans[:3])))
    return clusters, audit_rows


def cluster_anchors(anchors, window_bp):
    """Group anchors that sit close together on the SAME contig.

    Input:  every anchor from Phase 1 (both CONJscan machinery and integrases).
    Does:   sorts them per contig and starts a new cluster whenever the gap to
            the previous anchor is larger than `window_bp`. Because the anchors
            are CDS features, the resulting interval is already snapped to coding
            boundaries, which is what the spec asks for.
    Output: a list of lists of anchors, contig by contig, left to right.

    The biology behind clustering at all: the genes of a conjugative element are
    next to each other, in one block. Machinery genes scattered across a
    chromosome are unrelated leftovers, not one element.

    Caveat to know about: this is single-linkage clustering, so a chain of
    anchors each within `window_bp` of the next can grow much longer than the
    window itself. That is intentional (element genes are interrupted by cargo),
    and `--max-element-bp` is the guard against a chain that has run away; every
    such drop is written to the audit file with its measured span.
    """
    anchors_by_contig = {}
    for anchor in anchors:
        anchors_by_contig.setdefault(anchor["contig"], []).append(anchor)

    clusters = []
    for contig in sorted(anchors_by_contig):
        ordered = sorted(anchors_by_contig[contig], key=lambda a: (a["start"], a["end"]))
        current = [ordered[0]]
        current_max_end = ordered[0]["end"]

        for anchor in ordered[1:]:
            # Bases strictly between the running end of the cluster and this
            # anchor. Overlapping or touching anchors give a negative number or
            # zero, which is always within the window.
            gap = anchor["start"] - current_max_end - 1
            if gap <= window_bp:
                current.append(anchor)
                current_max_end = max(current_max_end, anchor["end"])
            else:
                clusters.append(current)
                current = [anchor]
                current_max_end = anchor["end"]

        clusters.append(current)

    return clusters


def merge_clusters_sharing_a_system(clusters, max_element_bp):
    """Re-join clusters that distance-clustering cut out of ONE CONJscan system.

    THE PROBLEM THIS FIXES. cluster_anchors groups anchors purely on how far
    apart they sit, so a conjugation system whose genes are spread a little wider
    than --window-bp gets reported as two or more separate elements. On the
    Phase 7 benchmark that was not an edge case: R391 was split into an 'ice' of
    56.6 kb plus a stray 'conjugative_region', because the gap between the two
    halves measured 15,010 bp against a 15,000 bp window - a TEN base pair
    margin on a threshold the spec itself calls convention rather than biology.
    Five of the eight genomes with calls showed the same signature: two rows,
    one sys_id.

    Two things go wrong when it happens. The element count is inflated, and -
    worse - the reported element is truncated, because `start`/`end` come from
    the surviving cluster alone. Merging R391's two halves takes the recovered
    span from 64% of the curated element to 94%.

    WHY MERGING IS SAFE RATHER THAN A LOOSER THRESHOLD. MacSyFinder has already
    decided these genes form one system, using its own co-localisation rules
    (counted in GENES between components, not in base pairs). Honouring that
    decision is not us relaxing a cutoff - it is us stopping overriding a
    judgement made with better information. Widening --window-bp instead would
    merge genuinely unrelated neighbours too.

    Input:  the distance-based clusters, each a list of anchors carrying sys_id.
    Does:   joins clusters that share any non-empty sys_id ON THE SAME CONTIG.
    Output: (merged clusters, audit rows describing every merge).

    Two deliberate restrictions:
      - anchors with an EMPTY sys_id (integrases found by Bakta product text,
        which belong to no CONJscan system) never cause a merge, or every
        unrelated integrase would pull the whole contig into one blob;
      - clusters on DIFFERENT contigs are never merged even when they share a
        sys_id. On a fragmented assembly MacSyFinder treats the proteome as one
        pseudo-replicon and can group hits across contigs; merging those would
        invent an interval spanning sequence that was never contiguous. That
        case is already reported through the spans_contigs flag.

    AND ONE SIZE GUARD, learned the hard way on the benchmark. A sys_id is not a
    promise of proximity: on the 10.1 Mb Streptomyces scabiei chromosome
    MacSyFinder put hits 2.19 Mb apart into one system, and merging them
    produced a 2,185,966 bp "element" that then blew past --max-element-bp and
    was DROPPED - turning a partial detection into a total miss. So a merge that
    would exceed max_element_bp is REFUSED and the original clusters are kept.
    The merge exists to repair an artefact of an arbitrary bp threshold, not to
    override the biological limit on how large an element can be.
    """
    # Group cluster indices by (contig, sys_id). A cluster carrying two system
    # IDs links both groups, which is why this needs a proper union rather than
    # a single pass.
    parent = list(range(len(clusters)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        root_i, root_j = find(i), find(j)
        if root_i != root_j:
            parent[max(root_i, root_j)] = min(root_i, root_j)

    first_cluster_with_key = {}
    for index, cluster in enumerate(clusters):
        contig = cluster[0]["contig"]
        for anchor in cluster:
            sys_id = (anchor.get("sys_id") or "").strip()
            if not sys_id:
                continue
            key = (contig, sys_id)
            if key in first_cluster_with_key:
                union(first_cluster_with_key[key], index)
            else:
                first_cluster_with_key[key] = index

    groups = {}
    for index in range(len(clusters)):
        groups.setdefault(find(index), []).append(index)

    merged = []
    audit_rows = []
    for root in sorted(groups):
        members = groups[root]
        joined = []
        for index in members:
            joined.extend(clusters[index])
        joined.sort(key=lambda a: (a["start"], a["end"]))

        # The size guard. Refuse a merge that would produce something larger than
        # any real element, and keep the pieces instead - a truncated call is far
        # more useful than an element dropped for being impossibly long.
        merged_span = max(a["end"] for a in joined) - min(a["start"] for a in joined) + 1
        if len(members) > 1 and merged_span > max_element_bp:
            shared = sorted({(a.get("sys_id") or "").strip()
                             for index in members for a in clusters[index]
                             if (a.get("sys_id") or "").strip()})
            audit_rows.append(audit_row(
                "NA",
                "evidence_recorded",
                "system_merge_refused_too_long",
                f"{len(members)} anchor clusters share the CONJscan system(s) "
                f"{', '.join(shared)} but merging them would span {merged_span} bp, "
                f"above the --max-element-bp threshold of {max_element_bp} bp. "
                "MacSyFinder can group hits that sit megabases apart on a large "
                "replicon; they are kept as separate clusters rather than joined "
                "into one implausible element.",
                contig=joined[0]["contig"],
                start=min(a["start"] for a in joined),
                end=max(a["end"] for a in joined),
            ))
            for index in sorted(members):
                merged.append(clusters[index])
            continue

        merged.append(joined)

        if len(members) > 1:
            shared = sorted({(a.get("sys_id") or "").strip()
                             for index in members for a in clusters[index]
                             if (a.get("sys_id") or "").strip()})
            spans = ", ".join(
                f"{min(a['start'] for a in clusters[i])}-{max(a['end'] for a in clusters[i])}"
                for i in sorted(members))
            audit_rows.append(audit_row(
                # `sample` is filled in by the caller, exactly as it is for the
                # integrase-attachment audit rows.
                "NA",
                "evidence_recorded",
                "clusters_merged_same_conjscan_system",
                f"{len(members)} anchor clusters ({spans}) were farther apart than "
                f"the clustering window but MacSyFinder assigned them to the same "
                f"system ({', '.join(shared)}), so they are reported as ONE element. "
                "Splitting them would both inflate the element count and truncate "
                "the element's reported span.",
                contig=joined[0]["contig"],
                start=min(a["start"] for a in joined),
                end=max(a["end"] for a in joined),
            ))

    return merged, audit_rows


def cluster_has_class(cluster, anchor_class):
    """True when at least one anchor in the cluster is of that class."""
    return any(anchor["anchor_class"] == anchor_class for anchor in cluster)


def cluster_looks_like_an_ime(cluster, has_mpf_system):
    """True when this cluster has the IME architecture: integrase + relaxase, no T4SS.

    An IME is mobilised by a helper element rather than moving itself, so by
    definition it encodes NO mating-pair apparatus. That is also why its
    machinery is SMALL - two genes, where an ICE carries a twenty-gene T4SS
    operon - and why the size floor below has to know which of the two it is
    looking at.

    `has_mpf_system` is passed in rather than derived here because the caller has
    already worked out whether CONJscan called a typed mating-pair system for
    this cluster's anchors; an IME must have neither the typed system nor an
    in-cluster mating-pair gene.
    """
    return (cluster_has_class(cluster, ANCHOR_INTEGRASE)
            and cluster_has_class(cluster, ANCHOR_RELAXASE)
            and not cluster_has_class(cluster, ANCHOR_T4SS)
            and not has_mpf_system)


def cluster_looks_like_an_aice(cluster, has_aice_system):
    """True when this cluster is an AICE: ICEscan called its AICE model here AND
    the anchors back it up with an integrase and its translocation machinery.

    BOTH halves are required, and the second is the important one. FtsK/SpoIIIE
    is a core chromosome-partitioning ATPase carried by essentially every
    bacterium - seeding an element on one would produce an "AICE" in every
    genome we ever run. Requiring ICEscan's assembled AICE model (translocase +
    a Rep protein + an integrase, within 10 genes of each other) is what makes
    the call specific; requiring the anchors too is what stops a model called
    somewhere else on the replicon from labelling this cluster.
    """
    return (has_aice_system
            and cluster_has_class(cluster, ANCHOR_INTEGRASE)
            and cluster_has_class(cluster, ANCHOR_AICE))


def keep_or_drop_cluster(cluster, min_element_bp, max_element_bp,
                         min_ime_element_bp=None, has_mpf_system=False,
                         has_aice_system=False):
    """Decide whether one cluster becomes a candidate element.

    Returns (keep, reason, detail). The three rejections, in the order they are
    tested:

    1. no relaxase and no T4SS component (spec §8 Phase 2, "keep clusters with
       >= relaxase or T4SS"). A cluster of integrases with no conjugation
       machinery is just a genome's normal complement of recombinases - every
       bacterium has several - and reporting each one as a mobile element would
       be pure noise.
    2. shorter than --min-element-bp. Real ICEs are tens of kilobases; a 3 kb
       span of machinery is a relic or a lone gene, not an element.
    3. longer than --max-element-bp. This is the runaway-chain guard described in
       cluster_anchors.

    The measured span is always reported in `detail`, because these thresholds
    are convention and the next reader may disagree with them.
    """
    start = min(anchor["start"] for anchor in cluster)
    end = max(anchor["end"] for anchor in cluster)
    span = end - start + 1

    has_relaxase = cluster_has_class(cluster, ANCHOR_RELAXASE)
    has_t4ss = cluster_has_class(cluster, ANCHOR_T4SS)
    is_aice = cluster_looks_like_an_aice(cluster, has_aice_system)

    # An AICE has NEITHER a relaxase NOR a mating-pair apparatus - that is what
    # it is, not what is wrong with it - so the conjugation test below would
    # throw away every one of them. It is admitted on its own evidence instead:
    # ICEscan's assembled AICE model plus the matching anchors.
    if not has_relaxase and not has_t4ss and not is_aice:
        classes = sorted({anchor["anchor_class"] for anchor in cluster})
        return False, "no_conjugation_anchor", (
            f"{len(cluster)} anchor(s) spanning {span} bp, of class(es) "
            f"{', '.join(classes)}: neither a relaxase nor a mating-pair (T4SS) "
            "component is present, so there is no evidence of conjugation and this "
            "is not seeded as an element."
        )

    # The size floor has to know WHICH KIND of element it is looking at.
    #
    # The floor is applied to the ANCHOR-CLUSTER SPAN, and that span scales with
    # the NUMBER of machinery genes, not with the element's length. An ICE
    # carries a ~20-gene T4SS operon, so its span is naturally tens of kb; an IME
    # carries a relaxase and an integrase - two genes - so its span is 1-6 kb.
    # An 8,000 bp floor therefore selects for ICEs by construction, which is
    # exactly what the Phase 7 IME pilot measured: six of twelve curated IMEs
    # were clustered and classified correctly and then dropped for size, Tn4451
    # at a span of 1,266 bp.
    #
    # So IME-architecture clusters get their own, lower floor. Note this is NOT
    # a general loosening: the ICE floor is untouched, and every one of the 18
    # ICE-pilot elements scores identically before and after.
    #
    # CAUTION FOR THE READER: the IME floor is an EMPIRICAL CUT, not a
    # biological constant. There is no published minimum size for an IME's
    # machinery that we could find, and the value here was chosen against a
    # handful of observations in one benchmark. Treat it as a knob, and say so
    # in any methods write-up rather than implying it is a property of IMEs.
    #
    # An AICE takes the same lower floor, for the same reason: its machinery is
    # three genes (translocase, Rep protein, integrase), not a twenty-gene
    # operon, so the ICE floor would select it out by construction.
    floor = min_element_bp
    floor_flag = "--min-element-bp"
    small_machinery = (cluster_looks_like_an_ime(cluster, has_mpf_system)
                       or is_aice)
    if min_ime_element_bp is not None and small_machinery:
        floor = min_ime_element_bp
        floor_flag = "--min-ime-element-bp"

    if span < floor:
        return False, "cluster_shorter_than_min", (
            f"machinery span is {span} bp, below the {floor_flag} threshold of "
            f"{floor} bp. Reported as too small to be a credible "
            "integrative element; the hits themselves are still in CONJscan's own "
            "output."
        )

    if span > max_element_bp:
        return False, "cluster_longer_than_max", (
            f"machinery span is {span} bp, above the --max-element-bp threshold of "
            f"{max_element_bp} bp. Anchors were chained together across the "
            "contig by the clustering window rather than forming one element."
        )

    return True, "", ""


# ── Confidence ───────────────────────────────────────────────────────────────

CONFIDENCE_RANK = {"low": 1, "medium": 2, "high": 3}


def assess_confidence(n_anchor_classes, machinery_intact, spans_contigs, at_contig_boundary,
                      boundary_method="none", require_trna_boundary=False,
                      from_profile_hits_only=False):
    """Give the candidate a confidence level, and say what lowered it.

    Returns (level, caps) where caps is a list of (level, reason, detail) - one
    per rule that fired, so the audit file can explain every "low" in the table.

    A call is only HIGH when all of these are true at once:
      * at least three of the four anchor classes are present, so it does not
        rest on a single hit;
      * every hit is on ONE contig - see below;
      * the machinery is complete (no low system wholeness, no decayed model, no
        truncated relaxase/VirB4);
      * the machinery does not run into a contig end, where the rest of the
        element would be invisible;
      * and, ONLY IF require_trna_boundary is on, the element's ends are known
        from a tRNA-anchored att pair.

    WHY THE BOUNDARY RULE IS OPTIONAL AND OFF BY DEFAULT
        The spec (§8 Phase 6) writes the rule as "high = 4 anchor classes +
        tRNA-anchored + single contig + intact", which folds two different
        questions into one number:

          1. how sure are we this IS an ICE?   (anchors, machinery, one contig)
          2. how sure are we WHERE IT ENDS?    (boundary_method)

        On a short-read assembly the second question usually has no answer - the
        flanks are simply not in the contig - so requiring a tRNA-anchored att
        pair makes `high` almost unreachable and drags down calls whose
        CLASSIFICATION is not in doubt. That is a poor trade: the module's main
        claim is "this AMR gene sits in a self-transmissible element", and that
        claim rests on question 1.

        So by default the boundary only informs, and the two facts are read side
        by side - confidence for the classification, boundary_method for the
        extent. Set require_trna_boundary to get the strict spec reading, which
        is the right choice on closed long-read assemblies where an unresolved
        boundary really is a warning sign rather than the norm.

        Note this is ONLY about the confidence label. Cargo assignment is
        unaffected either way: a de novo repeat never widens an element (see
        refine_candidate_boundaries), so an unresolved boundary always means the
        interval is the machinery span - a floor, never an invention.

    The contig rule is absolute, exactly as the spec requires: anything spanning
    contigs is capped at LOW no matter how good the rest of the evidence looks.
    MacSyFinder is run over the whole draft as one pseudo-replicon, so genes it
    joined across a contig break may simply be unrelated genes that happen to be
    adjacent in the file.
    """
    caps = []

    # The profile-hit fallback is the weakest evidence this module reports:
    # MacSyFinder saw the profiles but never assembled them into a system,
    # because some component of the quorum was absent. The row exists so a real
    # element is not silently lost, and it is capped at LOW so it can never be
    # mistaken for a system-backed call.
    if from_profile_hits_only:
        caps.append(("low", "no_assembled_system", (
            "MacSyFinder did not assemble these hits into a complete conjugation "
            "system, so the row rests on individual profile hits alone. See "
            "missing_components for which part of the machinery was not found; "
            "the element is reported because the remaining evidence was too "
            "strong to discard, not because a system was called."
        )))

    # boundary_method=None means NOT YET ASSESSED, and is what build_candidates
    # passes: it runs before Phase 3, so at that point no boundary has been looked
    # for. Judging it there would cap every element for a failure that has not
    # happened yet, then have to undo the cap - leaving a contradictory
    # "no att pair was found" line in the audit of an element whose att pair was
    # found moments later. finalise_confidence settles it once, afterwards.
    if require_trna_boundary and boundary_method is not None and boundary_method != "tRNA":
        if boundary_method == "denovo":
            caps.append(("medium", "boundary_denovo_only", (
                "the only candidate boundary is a de novo direct repeat, which on "
                "a real chromosome arises by chance often enough (~16% of "
                "arbitrary spans) that it was reported but not applied; the "
                "element's true extent is not established."
            )))
        else:
            caps.append(("medium", "boundary_not_resolved", (
                "no att pair was found, so the reported interval is the machinery "
                "span rather than the element's real ends - a floor, not a "
                "delimitation."
            )))

    if n_anchor_classes < MIN_ANCHOR_CLASSES_FOR_HIGH:
        caps.append(("medium", "few_anchor_classes", (
            f"only {n_anchor_classes} of the four anchor classes "
            f"(integrase / relaxase / coupling protein / T4SS) are present; "
            f"{MIN_ANCHOR_CLASSES_FOR_HIGH} are required for high confidence."
        )))

    if not machinery_intact:
        caps.append(("medium", "machinery_not_intact", (
            "the conjugation machinery is incomplete or degraded, so the element "
            "may no longer be able to do what its gene content suggests."
        )))

    if at_contig_boundary is None:
        caps.append(("medium", "contig_length_unknown", (
            "no length is known for this contig, so we cannot tell whether the "
            "element runs off the end of the assembly."
        )))
    elif at_contig_boundary:
        caps.append(("low", "at_contig_boundary", (
            "the machinery reaches a contig end: the element is very probably "
            "truncated by the assembly and its real extent is unknown."
        )))

    if spans_contigs:
        caps.append(("low", "spans_contigs", (
            "this system's hits are on more than one contig. MacSyFinder treats a "
            "draft assembly as a single ordered replicon, so genes on either side "
            "of a contig break can be joined into one system that does not exist. "
            "Capped at low confidence regardless of the other evidence."
        )))

    level = "high"
    for cap_level, _reason, _detail in caps:
        if CONFIDENCE_RANK[cap_level] < CONFIDENCE_RANK[level]:
            level = cap_level
    return level, caps


# ── Phase 4 assembled: clusters -> candidate rows ────────────────────────────

def build_candidates(sample, clusters, systems, contig_lengths,
                     min_element_bp, max_element_bp, boundary_bp,
                     min_ime_element_bp=None):
    """Turn the Phase 2 clusters into the deliverable rows.

    Input:  the clusters, the per-system summary from build_conjscan_anchors,
            and the contig lengths.
    Does:   applies the keep/drop rules, classifies what survives, works out the
            degradation, contig-boundary and contig-spanning flags, and sets a
            confidence level.
    Output: (rows, audit_rows). `rows` is what colocalise.py reads next.
    """
    rows = []
    audit_rows = []

    for cluster in clusters:
        contig = cluster[0]["contig"]
        start = min(anchor["start"] for anchor in cluster)
        end = max(anchor["end"] for anchor in cluster)

        # Which typed mating-pair systems this cluster's anchors belong to.
        # Needed BEFORE the keep/drop test, because the size floor depends on
        # whether this looks like an IME or an ICE.
        cluster_system_ids_early = {anchor["sys_id"] for anchor in cluster
                                    if anchor.get("sys_id")}
        has_mpf_system_early = any(
            systems.get(sys_id, {}).get("mpf_types", set())
            for sys_id in cluster_system_ids_early
            if systems.get(sys_id, {}).get("source") == SOURCE_CONJSCAN
        )
        has_aice_system_early = any(
            systems.get(sys_id, {}).get("icescan_class") == ICESCAN_MODEL_AICE
            for sys_id in cluster_system_ids_early
        )
        keep, reason, detail = keep_or_drop_cluster(
            cluster, min_element_bp, max_element_bp,
            min_ime_element_bp=min_ime_element_bp,
            has_mpf_system=has_mpf_system_early,
            has_aice_system=has_aice_system_early)
        if not keep:
            audit_rows.append(audit_row(
                sample, "dropped", reason, detail, contig=contig, start=start, end=end
            ))
            continue

        # --- which evidence is present -------------------------------------
        has_integrase = cluster_has_class(cluster, ANCHOR_INTEGRASE)
        has_relaxase = cluster_has_class(cluster, ANCHOR_RELAXASE)
        has_t4cp = cluster_has_class(cluster, ANCHOR_T4CP)
        has_t4ss = cluster_has_class(cluster, ANCHOR_T4SS)
        present_classes = [
            anchor_class for anchor_class in ANCHOR_CLASS_ORDER
            if cluster_has_class(cluster, anchor_class)
        ]

        # --- does the SYSTEM itself evidence a mating-pair apparatus? --------
        # `has_t4ss` above only says that SOME mating-pair profile was hit inside
        # this cluster. That alone must not make an ICE, because CONJscan's
        # relaxase-centred `MOB` model lists VirB4 as an ACCESSORY gene: a MOB
        # system is BY DEFINITION a relaxase-only system - it describes DNA that
        # can be picked up by someone else's machinery - and one incidental VirB4
        # inside it is not evidence that this cell can build a mating bridge.
        #
        # Only the typed models (`T4SS_type*`, and their decayed `dCONJ_type*`
        # counterparts) describe a complete mating-pair apparatus. So we ask, per
        # mating-pair hit, which model it was found under, and count only the
        # typed ones. Without this test a single accessory VirB4 would promote an
        # IME (mobilisable, needs a helper) straight to an ICE (predicted
        # self-transmissible) - the worst overcall this script could make, since
        # tier 6 is exactly the answer a regulator reads.
        # Ask the SYSTEM, not the individual hit. An earlier version of this test
        # looked at each mating-pair anchor's own model, which diverges from the
        # system view in a case that really happens: in T4SS_typeF both the
        # relaxase and the coupling protein are declared loner genes, so a typed
        # system can contribute those two while a separate MOB system in the same
        # cluster contributes the accessory VirB4. The per-hit test then said "no
        # apparatus" for a cluster that plainly had a typed T4SS system in it, and
        # wrote an audit line asserting something the row's own mpf_type column
        # contradicted. Whether CONJscan called a typed mating-pair SYSTEM here is
        # the question that matters, and mpf_types already answers it.
        # Which typed mating-pair systems do the anchors in this cluster belong
        # to? Computed here from the anchors' own system ids, because the tier
        # decision below needs the answer; the descriptive mpf_types list further
        # down is built the same way and reports it.
        cluster_system_ids = {anchor["sys_id"] for anchor in cluster
                              if anchor.get("sys_id")}

        # THE ONE PLACE THE TWO MODEL SETS ARE TREATED DIFFERENTLY.
        #
        # Every system-level property below - the mating-pair type, the
        # completeness, the decayed/truncated flags, and the set of contigs a
        # system's hits landed on - is read from CONJscan systems ONLY. ICEscan
        # contributes its hits' coordinates and its IME/AICE model class, and
        # nothing else. Two measured reasons:
        #
        #   * ICEscan's T4SS quorum is stricter and its model set is missing
        #     MOB.xml and the decayed dCONJ models, so its apparatus calls are
        #     not interchangeable with the ones every threshold here was tuned
        #     against. Keeping them out is what holds the tier-6 counts on the
        #     benchmark exactly where they were (self-transmissible 36,
        #     high-confidence 17, unchanged by adding ICEscan).
        #   * an ICEscan IME "system" is not a locus - both its mandatory genes
        #     are loners, so 17 of the 71 on this benchmark span more than 1,000
        #     genes. Its contig set would set spans_contigs on almost any draft
        #     assembly and cap perfectly good calls at low confidence.
        conjscan_system_ids = {
            sys_id for sys_id in cluster_system_ids
            if systems.get(sys_id, {}).get("source") == SOURCE_CONJSCAN
        }
        cluster_mpf_types = {
            mpf_type
            for sys_id in conjscan_system_ids
            for mpf_type in systems.get(sys_id, {}).get("mpf_types", set())
        }
        has_mpf_system = bool(cluster_mpf_types)

        # The classes ICEscan called for this cluster's anchors, if any.
        cluster_icescan_classes = {
            systems[sys_id]["icescan_class"]
            for sys_id in cluster_system_ids
            if systems.get(sys_id, {}).get("icescan_class")
        }

        # Does this cluster's machinery come from a system whose hits are spread
        # over more than one contig? Computed here rather than with the other
        # short-read flags further down, because the classification a few lines
        # below needs it. It is reported unchanged as the spans_contigs column.
        cluster_spans_contigs = any(
            len(systems.get(sys_id, {}).get("contigs", set())) > 1
            for sys_id in conjscan_system_ids
        )
        # ...but asking the system alone is not enough either, and the Phase 7
        # benchmark showed why. MacSyFinder LONER genes may sit anywhere on the
        # replicon, so a lone relaxase declared a loner of a typed T4SS model
        # drags that model's type letter across the whole chromosome. On
        # NC_013929 that produced a cluster holding exactly ONE integrase and ONE
        # relaxase - the textbook IME signature, has_t4cp=FALSE, has_t4ss=FALSE -
        # reported as `ice`, "predicted self-transmissible". The row contradicted
        # itself: missing_components read "coupling protein, mating-pair
        # apparatus" beside a tier-6 mobility claim, and the caller's own audit
        # had already refused to merge the real apparatus in, 1.2 Mb away.
        #
        # So a typed system is necessary but the apparatus must also BE here: at
        # least one mating-pair anchor in this cluster.
        #
        # THE EXCEPTION, and it is not hypothetical. On a fragmented assembly the
        # tra operon routinely lands on a different contig from the relaxase, so
        # requiring an in-cluster mating-pair gene would demote real ICEs the
        # moment the assembly breaks. Measured by rewriting the benchmark
        # genomes' coordinates into contigs while holding CONJscan's output
        # fixed: without this exemption, at 20 kb contigs two of the eighteen
        # ICE-pilot elements drop from `ice` to `ime` - one of them
        # ICE(Tn4371)6061, a spec section 8 Phase 7 named positive control. With
        # it, nothing changes under fragmentation and closed genomes are
        # unaffected. BacFlux's main input is short-read drafts, so this matters
        # more here than the closed-genome benchmark can show.
        has_mpf_apparatus = has_mpf_system and (has_t4ss or cluster_spans_contigs)

        if has_t4ss and not has_mpf_apparatus:
            audit_rows.append(audit_row(
                sample, "kept_flagged", "mpf_marker_without_typed_system",
                "a mating-pair marker was hit in this cluster, but CONJscan called "
                "no typed mating-pair system here - the hit came under a "
                "relaxase-centred MOB model, which lists VirB4 as an accessory "
                "gene. The marker is still reported "
                "in has_t4ss, but it was NOT counted as a mating-pair apparatus "
                "when classifying this cluster, so the element was not raised to "
                "ICE (predicted self-transmissible) on the strength of it. A "
                "relaxase-only system is mobilisable, not self-transmissible.",
                contig=contig, start=start, end=end,
            ))

        mge_class, mobility = classify_cluster(
            has_integrase, has_relaxase, has_mpf_apparatus
        )

        # --- what ICEscan's own model class adds -----------------------------
        # classify_cluster works purely from which anchors are present, and that
        # is still the rule for every conjugative class. ICEscan is allowed to
        # change the answer in exactly one situation - AICE - and is recorded
        # but not obeyed in the other.
        #
        # AICE: our anchor rules CANNOT reach this class, because they are built
        # around a relaxase and a mating bridge and an AICE has neither. There is
        # nothing for the ICEscan call to contradict, so it stands - and
        # cluster_looks_like_an_aice has already required the anchors to agree.
        #
        # IME: ICEscan's IME model is deliberately NOT allowed to promote
        # anything. Its quorum is two genes, either of which may be a loner
        # anywhere on the replicon, so an IME model firing near a cluster is weak
        # evidence about THIS locus. If our anchors say island (an integrase and
        # no relaxase) then there is no relaxase, and calling it an IME would
        # claim it can be mobilised by a helper - a real mobility claim, tier 5,
        # on evidence we do not have. Agreement and disagreement are both written
        # to the audit file so the comparison with ICEscan stays visible.
        mobility_tier_reason = TIER_ASSIGNED_DOWNSTREAM_REASON
        if ICESCAN_MODEL_AICE in cluster_icescan_classes and cluster_looks_like_an_aice(
                cluster, True):
            mge_class = MGE_CLASS_AICE
            mobility = AICE_MOBILITY
            mobility_tier_reason = AICE_TIER_REASON
            audit_rows.append(audit_row(
                sample, "kept_flagged", "class_taken_from_icescan_aice_model",
                "ICEscan assembled its AICE model here (FtsK/SpoIIIE translocase, "
                "a Rep protein and an integrase) and the anchors agree. Reported "
                "as an AICE: it moves as double-stranded DNA between hyphal "
                "compartments, so it is given NO mobility tier - tiers 5 and 6 "
                "are both conjugation and an AICE does neither. CAVEAT: the only "
                "AICE boundary we have validated (AICEScab56241) was recovered at "
                "MacSyFinder --coverage-profile 0.3, not the 0.5 used here; at 0.5 "
                "neither AICE call on the benchmark overlaps a curated one. Treat "
                "an AICE call as unvalidated.",
                contig=contig, start=start, end=end))
        elif ICESCAN_MODEL_IME in cluster_icescan_classes:
            agrees = mge_class == MGE_CLASS_IME
            audit_rows.append(audit_row(
                sample, "kept_flagged", "icescan_ime_model_agrees" if agrees
                else "icescan_ime_model_not_followed",
                f"ICEscan called its IME model over these anchors; our own rules "
                f"say '{mge_class}'."
                + (" The two agree." if agrees else
                   " OUR CALL STANDS. ICEscan's IME quorum is two genes, both "
                   "declared loners, so the model can fire from hits anywhere on "
                   "the replicon and says little about this locus. Following it "
                   "would assert the element is mobilisable by a helper - a "
                   "tier-5 mobility claim - which these anchors do not support."),
                contig=contig, start=start, end=end))

        # --- is the machinery actually intact? ------------------------------
        # Three independent ways of being broken, all reported by name so the
        # reader knows which one fired. Any of them downgrades the mobility
        # sentence itself, not just a flag column.
        # Reported in full (both tools) so a reader can trace every hit back,
        # but - as set out above - the health of the machinery is judged from
        # CONJscan's systems alone.
        contributing_systems = sorted({
            anchor["sys_id"] for anchor in cluster if anchor["sys_id"]
        })
        judged_systems = sorted(conjscan_system_ids)
        wholeness_values = [
            systems[sys_id]["wholeness"] for sys_id in judged_systems
            if systems.get(sys_id, {}).get("wholeness") is not None
        ]
        lowest_wholeness = min(wholeness_values) if wholeness_values else None

        degraded_reasons = []
        if lowest_wholeness is not None and lowest_wholeness < WHOLENESS_INTACT_MIN:
            degraded_reasons.append("low_system_wholeness")
        if any(systems.get(sys_id, {}).get("decayed") for sys_id in judged_systems):
            degraded_reasons.append("decayed_system_model")
        if any(systems.get(sys_id, {}).get("truncated") for sys_id in judged_systems):
            degraded_reasons.append("truncated_core_hit")
        machinery_intact = not degraded_reasons

        # An AICE is not a degraded ICE and must not be described as one: it has
        # no relaxase and no mating bridge BY DEFINITION, so "machinery
        # incomplete" would be an untrue sentence rather than a hedge.
        if mge_class == MGE_CLASS_AICE:
            machinery_intact = True
            degraded_reasons = []

        if not machinery_intact:
            mobility = mobility + DEGRADED_MOBILITY_SUFFIX
            audit_rows.append(audit_row(
                sample, "kept_flagged", "machinery_degraded",
                f"classified {mge_class}, but the machinery is degraded "
                f"({', '.join(degraded_reasons)}"
                + (f"; lowest sys_wholeness {lowest_wholeness}" if lowest_wholeness is not None else "")
                + f"). Mobility reported as '{mobility}'. Decayed elements are common "
                "and are where naive tools overcall.",
                contig=contig, start=start, end=end,
            ))

        # --- honest short-read flags ----------------------------------------
        spans_contigs = any(
            len(systems.get(sys_id, {}).get("contigs", set())) > 1
            for sys_id in judged_systems
        )

        contig_length = contig_lengths.get(contig)
        if contig_length is None:
            dist_to_start = None
            dist_to_end = None
            dist_to_nearest = None
            at_boundary = None
        else:
            dist_to_start = start - 1                  # bases of contig before it
            dist_to_end = max(contig_length - end, 0)  # bases of contig after it
            dist_to_nearest = min(dist_to_start, dist_to_end)
            at_boundary = dist_to_nearest <= boundary_bp

        # No boundary_method here on purpose - Phase 3 has not run yet, so the
        # element's ends are genuinely unknown at this point rather than
        # unresolved. finalise_confidence applies that rule once Phase 3 is done.
        confidence, caps = assess_confidence(
            len(present_classes), machinery_intact, spans_contigs, at_boundary,
            boundary_method=None,
            from_profile_hits_only=not contributing_systems,
        )
        for cap_level, cap_reason, cap_detail in caps:
            audit_rows.append(audit_row(
                sample, "kept_flagged", cap_reason,
                f"confidence capped at {cap_level}: {cap_detail}",
                contig=contig, start=start, end=end,
            ))

        # --- descriptive fields ---------------------------------------------
        # A single strand only means something when every anchor agrees; a real
        # element carries genes on both strands, so '.' is the usual answer and
        # downstream reads it as "unknown orientation".
        strands = {anchor["strand"] for anchor in cluster}
        strand = strands.pop() if len(strands) == 1 else "."

        relaxase_types = sorted({
            anchor["label"] for anchor in cluster
            if anchor["anchor_class"] == ANCHOR_RELAXASE
        })
        mpf_types = sorted({
            mpf_type
            for sys_id in judged_systems
            for mpf_type in systems.get(sys_id, {}).get("mpf_types", set())
        })
        # A T4SS_type*/dCONJ_type* model among the contributing systems means
        # CONJscan called a whole typed mating-pair system somewhere in this
        # cluster. Reported as a column so a reader can see at a glance whether
        # the call rests on a typed system; the classification decision itself
        # was already made above, per mating-pair hit, and audited there.
        mpf_typed_system = bool(mpf_types)

        integrase_products = [
            anchor["label"] for anchor in cluster
            if anchor["anchor_class"] == ANCHOR_INTEGRASE
        ]
        anchor_ids = [
            f"{anchor['feature_id']}({anchor['anchor_class']})" for anchor in cluster
        ]
        models = sorted({
            model
            for sys_id in contributing_systems
            for model in systems.get(sys_id, {}).get("models", set())
            if model
        })

        rows.append({
            "sample": sample,
            # Spec §9 identifier format: contig|type-start:end. Stable across
            # runs, and it is the join key a report can use.
            "mge_id": f"{contig}|{mge_class}-{start}:{end}",
            "mge_name": "NA",          # no curated naming database in this phase
            "element_type": ELEMENT_TYPE_FOR_CLASS[mge_class],
            "mge_class": mge_class,
            "mobility": mobility,
            # Always NA here: the tier belongs to an AMR gene, not to an element.
            # The reason column is what carries the meaning - and for an AICE it
            # says the tier is not merely unassigned but inapplicable.
            "mobility_tier": "NA",
            "mobility_tier_reason": mobility_tier_reason,
            "contig": contig,
            "start": str(start),
            "end": str(end),
            "length_bp": str(end - start + 1),
            "strand": strand,
            "n_anchors": str(len(cluster)),
            "n_anchor_classes": str(len(present_classes)),
            "anchor_classes": ",".join(present_classes),
            "has_integrase": _tsv_bool(has_integrase),
            "has_relaxase": _tsv_bool(has_relaxase),
            "has_t4cp": _tsv_bool(has_t4cp),
            "has_t4ss": _tsv_bool(has_t4ss),
            "relaxase_type": ",".join(relaxase_types) if relaxase_types else "NA",
            "mpf_type": ",".join(mpf_types) if mpf_types else "NA",
            "mpf_typed_system": _tsv_bool(mpf_typed_system),
            # ' | ' rather than ',' because product text is full of commas.
            "integrase_products": " | ".join(integrase_products) if integrase_products else "NA",
            "anchor_ids": ",".join(anchor_ids),
            "conjscan_systems": ",".join(contributing_systems) if contributing_systems else "NA",
            "conjscan_models": ",".join(models) if models else "NA",
            "evidence_sources": ",".join(sorted({
                anchor.get("source", SOURCE_CONJSCAN) for anchor in cluster
            })),
            "icescan_model_class": (",".join(sorted(cluster_icescan_classes))
                                    if cluster_icescan_classes else "NA"),
            "sys_wholeness_min": (
                f"{lowest_wholeness:g}" if lowest_wholeness is not None else "NA"
            ),
            "machinery_intact": _tsv_bool(machinery_intact),
            "degraded_reason": ",".join(degraded_reasons) if degraded_reasons else "NA",
            # A cluster whose anchors carry no sys_id came from the profile-hit
            # fallback: MacSyFinder found these profiles but never assembled them
            # into a system. Naming that in the row itself means a reader never
            # has to infer it from an empty conjscan_systems cell.
            "evidence_level": (
                "system" if contributing_systems else "profile_hits_only"
            ),
            # For an AICE, the relaxase and mating bridge are absent by
            # definition, so listing them as "missing" would read as a degraded
            # ICE. See AICE_MISSING_COMPONENTS.
            "missing_components": (
                AICE_MISSING_COMPONENTS if mge_class == MGE_CLASS_AICE
                else missing_components_for(cluster)
            ),
            # Phase 3 defaults. refine_candidate_boundaries() below overwrites
            # these - and start/end - when it finds a flanking att pair.
            "boundary_method": "none",
            "attL": "NA",
            "attR": "NA",
            "att_sequence": "NA",
            "att_length_bp": "0",
            "att_mismatches": "0",
            "att_trna": "NA",
            "machinery_start": str(start),
            "machinery_end": str(end),
            "contig_length": _text_or_na(contig_length),
            "dist_to_contig_start": _text_or_na(dist_to_start),
            "dist_to_contig_end": _text_or_na(dist_to_end),
            "dist_to_nearest_contig_end": _text_or_na(dist_to_nearest),
            "at_contig_boundary": _tsv_bool(at_boundary),
            "spans_contigs": _tsv_bool(spans_contigs),
            "confidence": confidence,
        })

    # Left-to-right along each contig: the same order a reader scans a genome in.
    rows.sort(key=lambda row: (row["contig"], int(row["start"]), int(row["end"])))
    return rows, audit_rows


def read_is_intervals(is_table_path):
    """Read the IS element table into {contig: [(start, end), ...]} for masking.

    Input:  {sample}_is_elements.tsv from isescan_table, or "" when the mobilome
            module ran without it.
    Output: dict of contig -> intervals. Empty dict when there is no table, which
            simply means the de novo att scan runs unmasked.

    Why the att search needs this at all: insertion sequences carry terminal
    repeats and duplicate a few bases of target DNA when they transpose, so an
    IS-rich neighbourhood is full of direct repeats that have nothing to do with
    ICE integration. Blanking them first is what stops those decoys outranking a
    real att site - the spec names it the most likely way to get Phase 3 wrong.
    """
    intervals_by_contig = {}
    if not is_table_path or not os.path.isfile(is_table_path):
        return intervals_by_contig

    with open(is_table_path, encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for record in reader:
            contig = (record.get("contig") or "").strip()
            start = att_search.to_int(record.get("start"))
            end = att_search.to_int(record.get("end"))
            if contig and start is not None and end is not None:
                intervals_by_contig.setdefault(contig, []).append((start, end))
    return intervals_by_contig


def refine_candidate_boundaries(sample, rows, genome_path, gff3_path,
                                is_intervals_by_contig, flank_window_bp,
                                min_element_bp, max_element_bp, boundary_bp):
    """Phase 3: replace each machinery span with the element's real ends.

    Input:  the candidate rows from build_candidates; the genome Bakta
            annotated; that GFF3 (for tRNAs); the IS intervals to mask, keyed by
            contig; and the search/size bounds.
    Does:   for every candidate, look for a flanking attL/attR pair
            (att_search.find_att_sites). When one is found, the row's start/end
            become the element those boundaries define and boundary_method
            records how - while machinery_start/machinery_end keep the original
            span, so nothing is lost and the two can always be compared.
    Output: (rows, audit_rows). Rows are modified in place and returned for
            convenience; every change, and every failure to find boundaries, is
            audited.

    Why replace start/end rather than only report the att sites: those columns
    are what colocalise.py intersects against the AMR genes to decide which are
    CARGO. An ICE is usually much larger than its tra cluster, so leaving the
    machinery span in place would systematically understate what travels with the
    element - the dangerous direction for an AMR report. The mge_id is
    deliberately NOT rewritten: it is a join key other tables may already
    reference, and silently repointing it would break them.

    Everything here degrades quietly. A missing genome, an unreadable GFF3 or a
    contig with no usable flanks simply leaves boundary_method='none' and the
    machinery span standing, which is the honest answer on a fragmented assembly
    where the element runs off the end of its contig.
    """
    audit_rows = []
    if not rows:
        return rows, audit_rows

    if not genome_path or not os.path.isfile(genome_path):
        audit_rows.append(audit_row(
            sample, "input_missing", "no_genome_for_att_search",
            "no genome FASTA was given, so element boundaries could not be "
            "resolved; the reported interval is the machinery span "
            "(boundary_method=none).",
        ))
        return rows, audit_rows

    sequences = att_search.read_fasta(genome_path)
    trnas_by_contig = att_search.group_by_contig(
        att_search.parse_trna_features(gff3_path)) if gff3_path else {}

    for row in rows:
        # An att site is the SCAR left by integrase-mediated site-specific
        # recombination: the element's attP and the host's attB recombine, and the
        # site is duplicated to attL and attR. No integrase means the element never
        # integrated, so there is no scar to find and any repeat the search turned
        # up would be something else entirely - a leftover IS end, a genuine
        # duplication, or noise. Searching anyway is not merely wasteful, it
        # manufactures boundaries that do not exist.
        #
        # This was caught on the KPNIH1 positive control: the two elements the
        # search "resolved" were both conjugative_region calls with no integrase
        # (one of them on a plasmid, which does not integrate at all), while the
        # one genuinely integrative element - the chromosomal IME - got nothing.
        # Exactly backwards, until this gate was added.
        if not _reads_true(row.get("has_integrase")):
            audit_rows.append(audit_row(
                sample, "not_applicable", "att_search_skipped_no_integrase",
                f"{row['mge_id']}: no integrase, so this element cannot have "
                "integrated and has no attL/attR scar to find. The reported "
                "interval is the machinery span.",
            ))
            continue

        contig = row["contig"]
        sequence = sequences.get(contig, "")
        if not sequence:
            audit_rows.append(audit_row(
                sample, "kept_flagged", "contig_not_in_genome_fasta",
                f"{row['mge_id']}: contig '{contig}' was not found in the genome "
                "FASTA, so no att search was possible.",
            ))
            continue

        machinery_start = int(row["machinery_start"])
        machinery_end = int(row["machinery_end"])
        result = att_search.find_att_sites(
            sequence, machinery_start, machinery_end,
            trnas=trnas_by_contig.get(contig, []),
            mask=is_intervals_by_contig.get(contig, []),
            flank_window_bp=flank_window_bp,
            min_element_bp=min_element_bp,
            max_element_bp=max_element_bp,
        )

        row["boundary_method"] = result["boundary_method"]
        row["attL"] = result["att_left"]
        row["attR"] = result["att_right"]
        row["att_sequence"] = result["att_sequence"]
        row["att_length_bp"] = str(result["att_length_bp"])
        row["att_mismatches"] = str(result["att_mismatches"])
        row["att_trna"] = result["trna"]

        if result["boundary_method"] == "none":
            audit_rows.append(audit_row(
                sample, "kept_flagged", "element_boundaries_not_resolved",
                f"{row['mge_id']}: no flanking att pair was found, so the "
                f"reported interval is the machinery span "
                f"({machinery_start}-{machinery_end}). The real element is "
                "probably larger; on a fragmented assembly the flanks are often "
                "simply not present.",
            ))
            continue

        # A DE NOVO boundary is REPORTED but never ACTED ON. This is the single
        # most important restraint in Phase 3, so here is the measurement behind
        # it, taken on the KPNIH1 chromosome (5.4 Mb) with the real IS mask:
        #
        #   300 randomly placed 15 kb spans, none of them an ICE
        #   -> 22% came back with a confident de novo "boundary"
        #   -> 16% still did after the repeat-family guard was added
        #
        # A one-in-six error rate is perfectly acceptable for a HYPOTHESIS and
        # completely unacceptable for something that silently redefines the
        # element. When start/end are widened, colocalise.py treats every AMR gene
        # in the new interval as CARGO of the element, so a fabricated 50 kb
        # boundary turns unrelated chromosomal genes into "predicted
        # self-transmissible" - and the genes most often swallowed are exactly the
        # ones that must never be called mobile: a gyrA point mutation is the
        # textbook INTRINSIC determinant, and it sits on the chromosome where
        # these spurious repeats live.
        #
        # Mode A (tRNA-anchored) is a different proposition: it starts from a
        # position integrases are known to target and requires the probe to be the
        # 3' end of an actual annotated tRNA, so it is specific enough to act on.
        #
        # The de novo columns stay in the output because they are a real lead for
        # a human to follow - which is what the spec's boundary_method column is
        # for. They just do not move the element.
        if result["boundary_method"] == "denovo":
            audit_rows.append(audit_row(
                sample, "kept_flagged", "denovo_att_reported_not_applied",
                f"{row['mge_id']}: a de novo direct repeat "
                f"({result['att_length_bp']} bp) was found at "
                f"{result['att_left']} / {result['att_right']}, implying an "
                f"element of {result['element_length_bp']} bp. It is REPORTED but "
                f"NOT used as the element boundary: on a real chromosome ~16% of "
                f"arbitrary spans yield such a repeat by chance, so the interval "
                f"stays the machinery span ({machinery_start}-{machinery_end}) "
                "and cargo is not assigned from it. Treat the att columns as a "
                "lead to check by hand, not as a delimitation.",
            ))
            continue

        # Boundaries found by the tRNA-anchored mode: widen the element to what
        # actually travels.
        row["start"] = str(result["element_start"])
        row["end"] = str(result["element_end"])
        row["length_bp"] = str(result["element_length_bp"])
        added_bp = result["element_length_bp"] - (machinery_end - machinery_start + 1)

        # Everything derived from the interval must be recomputed, or the row
        # ends up describing two different elements at once. build_candidates
        # worked these out from the MACHINERY span; the element is now wider, so
        # it may reach a contig end the machinery did not come near. Leaving them
        # stale would be quietly dangerous rather than merely untidy:
        # at_contig_boundary feeds assess_confidence, so a widened element
        # running off the end of its contig would keep a high confidence it no
        # longer deserves - the report claiming a complete element on evidence
        # that has just been truncated by the assembly.
        #
        # The CONFIDENCE that follows from these flags is settled later, in one
        # pass, by finalise_confidence - see the note there.
        contig_length = att_search.to_int(row.get("contig_length"))
        if contig_length is not None:
            dist_to_start = result["element_start"] - 1
            dist_to_end = max(contig_length - result["element_end"], 0)
            dist_to_nearest = min(dist_to_start, dist_to_end)
            at_boundary = dist_to_nearest <= boundary_bp
            row["dist_to_contig_start"] = str(dist_to_start)
            row["dist_to_contig_end"] = str(dist_to_end)
            row["dist_to_nearest_contig_end"] = str(dist_to_nearest)
            row["at_contig_boundary"] = _tsv_bool(at_boundary)
        audit_rows.append(audit_row(
            sample, "boundaries_resolved", f"att_pair_found_{result['boundary_method']}",
            f"{row['mge_id']}: attL {result['att_left']} / attR "
            f"{result['att_right']} ({result['att_length_bp']} bp repeat, "
            f"{result['att_mismatches']} mismatch(es)"
            + (f", at {result['trna']}" if result["trna"] not in ("", "NA") else "")
            + f"). The element is {result['element_length_bp']} bp, "
            f"{added_bp} bp more than the machinery span alone; the extra is "
            "cargo that travels with it.",
        ))

    return rows, audit_rows


def finalise_confidence(sample, rows, require_trna_boundary=False):
    """Settle every candidate's confidence ONCE, after Phase 3 has run.

    WHY THIS IS A SEPARATE PASS
        Confidence depends on how well the element's ENDS are known, and that is
        only decided in Phase 3 - after build_candidates has already produced the
        rows. Two orderings were possible and one of them is wrong:

          * judge the boundary inside build_candidates, then correct it later.
            That caps every element for a failure that has not happened yet and
            writes "no att pair was found" into the audit of elements whose att
            pair is found seconds later. The audit is the thing a reader trusts
            to explain the table, so a contradictory line in it is worse than no
            line at all.
          * leave the boundary out of the first assessment (boundary_method=None)
            and settle it here, once, when the answer is actually known.

        This is the second one.

    Takes in: the candidate rows, already widened (or not) by Phase 3, with
              boundary_method and the contig-distance flags final.
    Does:     re-runs assess_confidence with the real boundary_method and writes
              one audit line per cap, so every non-high call has a stated reason.
    Output:   (rows, audit_rows); rows are modified in place.
    """
    audit_rows = []
    for row in rows:
        confidence, caps = assess_confidence(
            att_search.to_int(row.get("n_anchor_classes")) or 0,
            _reads_true(row.get("machinery_intact")),
            _reads_true(row.get("spans_contigs")),
            # Tri-state on purpose: NA here means "no contig length was known",
            # which caps confidence, and is NOT the same as FALSE.
            _reads_tristate(row.get("at_contig_boundary")),
            boundary_method=row.get("boundary_method") or "none",
            require_trna_boundary=require_trna_boundary,
            # Re-read from the row: this pass rebuilds the confidence from
            # scratch, so a cap that is not restated here is silently dropped.
            # That is exactly how the profile-hit fallback first came out as
            # 'high' - the cap was applied in build_candidates and then thrown
            # away the moment boundaries were settled.
            from_profile_hits_only=(
                row.get("evidence_level") == "profile_hits_only"
            ),
        )
        if confidence != row["confidence"]:
            audit_rows.append(audit_row(
                sample, "kept_flagged", "confidence_settled_after_boundary_search",
                f"{row['mge_id']}: confidence moves from {row['confidence']} to "
                f"{confidence} now that the boundary search has run "
                f"(boundary_method={row.get('boundary_method') or 'none'})"
                + (": " + "; ".join(detail for _lvl, _reason, detail in caps)
                   if caps else ""),
                contig=row.get("contig"), start=row.get("start"), end=row.get("end"),
            ))
            row["confidence"] = confidence
    return rows, audit_rows


# ── The command-line interface ───────────────────────────────────────────────

def build_parser():
    """Define the CLI. Every threshold is exposed, because all three are
    conventions rather than biology and a user working on, say, Bacteroidetes ICEs
    may reasonably want different ones."""
    parser = argparse.ArgumentParser(
        description=(
            "Turn CONJscan machinery hits plus Bakta annotation into ICE/IME "
            "candidate rows for the mobilome co-localisation step. Implements "
            "spec section 8 phases 0, 1, 2, 3, 4 and 6. The att-site search "
            "(phase 3) runs in every mode - the constraint is assembly "
            "contiguity, not the sequencer - and degrades to "
            "boundary_method=none when the flanks are absent. Only a "
            "tRNA-anchored boundary is applied to the element interval; a de "
            "novo repeat is reported but not acted on."
        )
    )
    parser.add_argument("--sample", required=True,
                        help="Sample name, written into every output row.")
    parser.add_argument("--conjscan-tsv", default=None,
                        help="CONJscan/MacSyFinder best_solution.tsv, or the "
                             "MacSyFinder output directory containing it. May be "
                             "absent: most isolates have no conjugative system.")
    parser.add_argument("--min-ime-element-bp", type=int,
                        default=DEFAULT_MIN_IME_ELEMENT_BP,
                        help="Size floor for clusters with the IME architecture "
                             "(integrase + relaxase, no mating-pair apparatus), "
                             "whose machinery is two genes rather than a "
                             "twenty-gene operon. ICE-architecture clusters keep "
                             f"--min-element-bp. Default {DEFAULT_MIN_IME_ELEMENT_BP}. "
                             "An empirical cut, not a published bound.")
    parser.add_argument("--icescan-tsv", default=None,
                        help="OPTIONAL second MacSyFinder best_solution.tsv, run "
                             "with the ICEscan model set (or its output "
                             "directory). ICEscan is a fork of CONJScan 2.0.1 by "
                             "the same authors: it adds an IME model, an AICE "
                             "model and 21 profiles, but drops MOB, the decayed "
                             "system models and the whole Plasmids set, so it is "
                             "unioned with CONJscan rather than replacing it. "
                             "Only its integrase anchors, its Gram-positive "
                             "relaxase families and its IME/AICE model classes "
                             "are used. ABSENT BY DEFAULT: without it the result "
                             "is exactly what CONJscan alone produces. Its models "
                             "are CC BY-NC-SA, so they are fetched at runtime and "
                             "never shipped with BacFlux.")
    parser.add_argument("--icescan-hmmer-dir", default=None,
                        help="ICEscan's hmmer_results/ directory, used the same "
                             "way --conjscan-hmmer-dir is: only when no system "
                             "was assembled at all. Defaults to hmmer_results/ "
                             "beside the ICEscan best_solution.tsv.")
    parser.add_argument("--conjscan-hmmer-dir", default=None,
                        help="MacSyFinder's hmmer_results/ directory. Read ONLY "
                             "when no complete system was assembled, to recover an "
                             "element that failed the quorum by one component "
                             "(reported at low confidence, with the missing part "
                             "named). Defaults to hmmer_results/ beside "
                             "best_solution.tsv.")
    parser.add_argument("--bakta-gff", required=True,
                        help="Bakta GFF3 for the same sample. Supplies the "
                             "coordinates of each CONJscan protein hit and the "
                             "integrase anchors.")
    parser.add_argument("--contig-lengths", default=None,
                        help="Optional contig<TAB>length table (or a .fai). When "
                             "omitted, lengths are read from the GFF3's own "
                             "##sequence-region lines.")
    parser.add_argument("--window-bp", type=int, default=DEFAULT_WINDOW_BP,
                        help=f"Largest gap between two anchors that still counts as "
                             f"the same element (default {DEFAULT_WINDOW_BP}).")
    parser.add_argument("--min-element-bp", type=int, default=DEFAULT_MIN_ELEMENT_BP,
                        help=f"Candidates spanning less than this are dropped "
                             f"(default {DEFAULT_MIN_ELEMENT_BP}).")
    parser.add_argument("--max-element-bp", type=int, default=DEFAULT_MAX_ELEMENT_BP,
                        help=f"Candidates spanning more than this are dropped as "
                             f"runaway anchor chains (default {DEFAULT_MAX_ELEMENT_BP}).")
    parser.add_argument("--boundary-bp", type=int, default=DEFAULT_BOUNDARY_BP,
                        help=f"A candidate within this many bases of a contig end "
                             f"is flagged as probably truncated and capped at low "
                             f"confidence (default {DEFAULT_BOUNDARY_BP}).")
    parser.add_argument("--genome", default="",
                        help="genome FASTA Bakta annotated; enables the Phase 3 "
                             "att-site search that resolves element boundaries")
    parser.add_argument("--is-table", default="",
                        help="IS element table, masked out before the de novo att "
                             "scan so transposon repeats cannot masquerade as att sites")
    parser.add_argument("--att-flank-window-bp", type=int,
                        default=att_search.DEFAULT_FLANK_WINDOW_BP,
                        help="how far beyond the machinery span to look for the "
                             "element's real ends")
    parser.add_argument("--require-trna-boundary-for-high", action="store_true",
                        help="Apply the spec's strict §8 Phase 6 rule: an element "
                             "may only reach HIGH confidence if a tRNA-anchored "
                             "att pair fixed its ends. OFF by default, because on "
                             "a fragmented short-read assembly the flanks are "
                             "usually absent, so the rule would cap almost every "
                             "call regardless of how good the machinery evidence "
                             "is - it conflates 'is this an ICE?' with 'where "
                             "does it end?'. Worth turning ON for closed "
                             "long-read assemblies, where an unresolved boundary "
                             "is genuinely a warning rather than the norm. Either "
                             "way boundary_method is reported, and cargo is never "
                             "assigned from an unresolved boundary.")
    parser.add_argument("--replicons", default="",
                        help="Per-contig replicon table from rule "
                             "mobilome_replicons. Used for one definitional "
                             "check: an ICE integrates into a CHROMOSOME, so an "
                             "element with full machinery on a PLASMID contig is "
                             "a conjugative plasmid region, not an ICE. Optional "
                             "- without it that check is skipped and the call is "
                             "left as it was.")
    parser.add_argument("--integrase-window-bp", type=int,
                        default=DEFAULT_INTEGRASE_WINDOW_BP,
                        help="How far from a machinery cluster an integrase may "
                             "sit and still belong to the same element. Much "
                             "larger than --window-bp on purpose: machinery is an "
                             "operon and clusters tightly, while the integrase "
                             "sits at the element boundary. Default "
                             f"{DEFAULT_INTEGRASE_WINDOW_BP}.")
    parser.add_argument("--out-table", required=True,
                        help="Candidate element TSV (read by colocalise.py).")
    parser.add_argument("--out-audit", required=True,
                        help="Decision trail TSV: one row per drop or downgrade.")
    return parser


def summarise(sample, rows, audit_rows, out_audit):
    """Build the one-line summary printed to the run log.

    Not used for filtering by anything; it just makes the result visible without
    opening a file. It always names the classes separately, because "3 elements
    found" would be a misleading thing to read for a sample whose three elements
    are all passive islands.
    """
    counts = {
        MGE_CLASS_ICE: 0,
        MGE_CLASS_IME: 0,
        MGE_CLASS_ISLAND: 0,
        MGE_CLASS_CONJ_REGION: 0,
        MGE_CLASS_AICE: 0,
    }
    for row in rows:
        counts[row["mge_class"]] += 1

    confidences = {"high": 0, "medium": 0, "low": 0}
    for row in rows:
        confidences[row["confidence"]] += 1

    n_dropped = sum(1 for row in audit_rows if row["action"] == "dropped")

    return (
        f"Sample {sample}: {len(rows)} candidate element(s) - "
        f"{counts[MGE_CLASS_ICE]} ICE (predicted self-transmissible), "
        f"{counts[MGE_CLASS_IME]} IME (mobilisable with a helper), "
        f"{counts[MGE_CLASS_ISLAND]} passive island, "
        f"{counts[MGE_CLASS_CONJ_REGION]} unbounded conjugative region, "
        f"{counts[MGE_CLASS_AICE]} AICE (no conjugation tier); "
        f"confidence high {confidences['high']}, medium {confidences['medium']}, "
        f"low {confidences['low']}; {n_dropped} cluster(s) dropped "
        f"(reasons in {out_audit})."
    )


def main(argv=None):
    """Run the whole thing: read, anchor, cluster, classify, write.

    Exit codes: 0 whenever the biology says "nothing to report" (no CONJscan
    file, no systems, no annotation) - those write a well-formed empty table and
    an audit line saying why. 1 only for a wiring or format problem, where
    continuing would mean publishing a wrong answer.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    audit_rows = []
    rows = []

    def finish(return_code=0):
        """Write both outputs and print the summary. Called from every exit
        path, so an empty result is still a complete, readable pair of files."""
        write_tsv(args.out_table, OUTPUT_COLUMNS, rows)
        write_tsv(args.out_audit, AUDIT_COLUMNS, audit_rows)
        print(summarise(args.sample, rows, audit_rows, args.out_audit))
        return return_code

    # --- the annotation, without which nothing can be placed ----------------
    if not args.bakta_gff or not os.path.exists(args.bakta_gff):
        _warn(
            f"Bakta GFF3 not found at '{args.bakta_gff}'. Without it a CONJscan "
            "protein hit cannot be given coordinates, so no element can be called. "
            "Writing empty (but complete) outputs."
        )
        audit_rows.append(audit_row(
            args.sample, "input_missing", "bakta_gff_missing",
            f"no Bakta GFF3 at '{args.bakta_gff}'. CONJscan reports protein IDs "
            "only, so without the annotation there are no coordinates to cluster "
            "and no integrase anchors. No ICE/IME call was attempted for this "
            "sample - this is 'not looked at', not 'looked at and found nothing'.",
        ))
        return finish(0)

    features_by_id, cds_features, gff_contig_lengths = parse_bakta_gff(args.bakta_gff)
    if not cds_features:
        _warn(f"no CDS features found in '{args.bakta_gff}'.")
        audit_rows.append(audit_row(
            args.sample, "input_missing", "bakta_gff_has_no_cds",
            f"'{args.bakta_gff}' parsed but contained no CDS rows, so there is "
            "nothing to place hits against and no products to search for "
            "integrases.",
        ))
        return finish(0)

    # --- contig lengths: the supplied table wins, the GFF3 is the fallback ---
    contig_lengths = read_contig_lengths(args.contig_lengths)
    if not contig_lengths:
        contig_lengths = gff_contig_lengths
        if args.contig_lengths:
            audit_rows.append(audit_row(
                args.sample, "input_missing", "contig_lengths_unusable",
                f"'{args.contig_lengths}' was given but no contig/length pairs could "
                "be read from it; fell back to the GFF3's ##sequence-region lines.",
            ))

    # --- the CONJscan result ------------------------------------------------
    conjscan_path = resolve_conjscan_path(args.conjscan_tsv)
    if conjscan_path is None:
        audit_rows.append(audit_row(
            args.sample, "input_missing", "conjscan_output_missing",
            f"no CONJscan result at '{args.conjscan_tsv}'. Either CONJscan was not "
            "run (it is opt-in: its models are CC BY-NC-SA licensed) or it wrote "
            "nothing. No conjugation machinery was assessed, so any AMR gene on the "
            "chromosome stays at the lowest mobility tier for want of evidence, "
            "not because the evidence was negative.",
        ))
        return finish(0)

    try:
        hits = read_conjscan_hits(conjscan_path)
    except ValueError as error:
        # A header we do not recognise is a tool-version change. Stopping is the
        # right answer: silently reporting "no conjugation machinery" from a file
        # we failed to parse would be a wrong result that looks like a real one.
        sys.stderr.write(f"ERROR: {error}\n")
        return 1

    # --- the optional second model set --------------------------------------
    # ICEscan is unioned in here, at the hit level, because both tools are
    # MacSyFinder and write the same 22 columns - so everything downstream sees
    # one hit list and does not need to know there were two runs. What each tool
    # is allowed to CONCLUDE is decided later, per rule, from the `source` tag
    # that read_conjscan_hits put on every hit.
    #
    # When --icescan-tsv is not given, `icescan_hits` stays empty and this whole
    # block is a no-op: the caller behaves exactly as it did before ICEscan
    # existed. That equivalence is deliberate and is pinned by a test.
    icescan_hits = []
    icescan_path = resolve_conjscan_path(args.icescan_tsv)
    if args.icescan_tsv and icescan_path is None:
        audit_rows.append(audit_row(
            args.sample, "input_missing", "icescan_output_missing",
            f"--icescan-tsv was given as '{args.icescan_tsv}' but no "
            "best_solution.tsv could be read there. Continuing with CONJscan "
            "alone: no IME or AICE model class and no ICEscan integrase anchors "
            "were available, so those elements can only be found if CONJscan's "
            "own models reach them.",
        ))
    elif icescan_path is not None:
        try:
            icescan_hits = read_conjscan_hits(icescan_path, source=SOURCE_ICESCAN)
        except ValueError as error:
            sys.stderr.write(f"ERROR: {error}\n")
            return 1
        if not icescan_hits:
            icescan_hits = read_conjscan_profile_hits(
                args.icescan_hmmer_dir or default_hmmer_dir(icescan_path),
                source=SOURCE_ICESCAN)
        if icescan_hits:
            audit_rows.append(audit_row(
                args.sample, "evidence_recorded", "icescan_hits_unioned",
                f"{len(icescan_hits)} hit(s) from the ICEscan model set were "
                f"unioned with CONJscan's {len(hits)}. ICEscan contributes "
                "integrase anchors, the Gram-positive and IME relaxase families "
                "CONJScan 2.1.0 does not model, and its IME/AICE model classes. "
                "It does NOT contribute mating-pair apparatus calls or element "
                "spans - its T4SS quorum is stricter and its IME systems are not "
                "loci - so tier-6 (predicted self-transmissible) calls rest on "
                "CONJscan alone.",
            ))
    # CONJscan's own quorum fallback, keyed on CONJSCAN'S result alone. It has to
    # be asked before the union, not after: if ICEscan happened to report one
    # unrelated hit, an `if not hits` test over the merged list would silently
    # skip this recovery and lose the element it exists for.
    if not hits:
        # No SYSTEM was assembled. Before accepting that as "no machinery", look
        # at the individual profile hits: MacSyFinder needs a full quorum, and a
        # genuine element missing one component of it would otherwise disappear
        # entirely (ICEVflInd1 on the Phase 7 benchmark - twelve mating-pair
        # profiles, a coupling protein and VirB4, but no relaxase, so nothing was
        # reported at all). Anything found this way is labelled
        # evidence_level=profile_hits_only and capped at low confidence.
        hits = read_conjscan_profile_hits(args.conjscan_hmmer_dir or
                                          default_hmmer_dir(conjscan_path))
        if hits:
            profiles = sorted({hit["gene_name"] for hit in hits})
            audit_rows.append(audit_row(
                args.sample, "evidence_recorded", "no_system_using_profile_hits",
                f"CONJscan assembled no complete system from '{conjscan_path}', but "
                f"{len(hits)} individual profile hit(s) passed its own e-value and "
                f"coverage thresholds ({', '.join(profiles)}). These are reported as "
                "a weaker evidence tier - evidence_level=profile_hits_only, "
                "confidence capped at low, and missing_components naming the part "
                "of the machinery that was absent - rather than discarded, because "
                "a real element can fail the quorum by one component.",
            ))
        elif not icescan_hits:
            # Nothing from either model set. Only now is "no machinery" the
            # honest answer.
            audit_rows.append(audit_row(
                args.sample, "input_missing", "conjscan_found_no_systems",
                f"'{conjscan_path}' contains no system rows and no profile hits "
                "were found either: CONJscan ran and found no conjugation "
                "machinery. This is the expected result for most environmental "
                "isolates and is a genuine negative, not a failure.",
            ))
            return finish(0)

    # The union everything downstream works from. Done AFTER the fallback above
    # so that CONJscan's recovery is judged on CONJscan's evidence, and BEFORE
    # the anchor phase so that an element only ICEscan can see - a Gram-positive
    # IME whose relaxase family CONJScan 2.1.0 does not model - is not lost to an
    # early return.
    hits = hits + icescan_hits

    # --- Phase 1: anchors ---------------------------------------------------
    machinery_anchors, hmm_integrases, systems, machinery_audit = build_conjscan_anchors(
        args.sample, hits, features_by_id
    )
    audit_rows.extend(machinery_audit)

    integrase_anchors, integrase_audit = find_integrase_anchors(args.sample, cds_features)
    audit_rows.extend(integrase_audit)

    # One pool of integrase candidates from both sources, product-text first so
    # that it is the one kept when the same CDS is found twice.
    #
    # THE SAME CDS OFTEN IS found twice: an ICEscan Phage_integrase hit and a
    # Bakta "tyrosine recombinase XerC" product line can be the same gene. Two
    # anchors on one feature would count it twice in n_anchors and list it twice
    # in anchor_ids. The product-text anchor is the one kept, because its label
    # is a sentence a reader can check rather than a profile accession - and
    # the corroboration is not lost, it is recorded in the audit file.
    corroborated = 0
    product_text_by_feature = {anchor["feature_id"]: anchor
                               for anchor in integrase_anchors}
    for anchor in hmm_integrases:
        twin = product_text_by_feature.get(anchor["feature_id"])
        if twin is None:
            integrase_anchors.append(anchor)
            continue
        corroborated += 1
        # Keep the product-text anchor, but do not throw away what ICEscan knew
        # about it. A product-text anchor carries no system, so the ICEscan
        # system id and model come across - that is how the IME/AICE model class
        # reaches build_candidates for an integrase both sources found. Dropping
        # it would mean the BETTER-evidenced integrase (two independent methods
        # agreeing) silently lost the model class a worse-evidenced one keeps.
        if not twin["sys_id"]:
            twin["sys_id"] = anchor["sys_id"]
            twin["model_fqn"] = anchor["model_fqn"]
    if corroborated:
        audit_rows.append(audit_row(
            args.sample, "evidence_recorded", "integrase_corroborated_by_icescan",
            f"{corroborated} integrase(s) were found by BOTH the Bakta product "
            "text and an ICEscan integrase profile. Counted once, keeping the "
            "Bakta product text as the label because it names the protein in "
            "words. Agreement between an annotation and an HMM is the strongest "
            "integrase evidence this module has.",
        ))

    if not machinery_anchors:
        audit_rows.append(audit_row(
            args.sample, "input_missing", "no_machinery_hit_could_be_placed",
            f"CONJscan reported {len(hits)} hit(s) but none of their protein IDs "
            "were found in the Bakta GFF3, so nothing could be given coordinates. "
            "Check that CONJscan was run on this sample's Bakta .faa.",
        ))
        return finish(0)

    # --- Phase 2 + Phase 4 --------------------------------------------------
    # Cluster the MACHINERY only, then attach integrases with a wider radius -
    # see attach_nearby_integrases for why the two anchor types need different
    # windows, and what it cost us when they shared one.
    clusters = cluster_anchors(machinery_anchors, args.window_bp)
    # Undo any split that distance-clustering made THROUGH a single CONJscan
    # system. Runs before integrases are attached, because integrases carry no
    # sys_id and so must not influence which clusters belong together.
    clusters, merge_audit = merge_clusters_sharing_a_system(
        clusters, args.max_element_bp)
    for row in merge_audit:
        row["sample"] = args.sample
    audit_rows.extend(merge_audit)
    # The att probe for the integrase tie-break (FIX 4). When two or more
    # integrases are in the running for one cluster - which is the normal case
    # once ICEscan's hits join the pool - the one that actually produces a
    # flanking attL/attR pair is the one that put this element here. This closure
    # answers "would including an integrase at this span yield an att pair?" by
    # running the same Phase 3 search that will later fix the boundaries.
    #
    # Only a tRNA-ANCHORED pair counts, exactly as in refine_candidate_boundaries:
    # a de novo repeat turns up in ~16% of arbitrary chromosomal spans by chance,
    # which is far too noisy to decide anything with.
    #
    # Without a genome FASTA there is nothing to search, so the probe is None and
    # the ranking falls back to distance then source.
    # The IS intervals are read once here and reused by Phase 3 below: the att
    # search must mask insertion sequences out, or their terminal inverted
    # repeats flood the candidate att pairs.
    is_intervals_by_contig = read_is_intervals(args.is_table)

    att_probe = None
    if args.genome and os.path.isfile(args.genome):
        probe_sequences = att_search.read_fasta(args.genome)
        probe_trnas = att_search.group_by_contig(
            att_search.parse_trna_features(args.bakta_gff))

        def att_probe(contig, span_start, span_end):
            sequence = probe_sequences.get(contig, "")
            if not sequence:
                return False
            result = att_search.find_att_sites(
                sequence, span_start, span_end,
                trnas=probe_trnas.get(contig, []),
                mask=is_intervals_by_contig.get(contig, []),
                flank_window_bp=args.att_flank_window_bp,
                min_element_bp=args.min_element_bp,
                max_element_bp=args.max_element_bp,
            )
            return result["boundary_method"] == "tRNA"

    clusters, integrase_attach_audit = attach_nearby_integrases(
        clusters, integrase_anchors, args.integrase_window_bp,
        att_supports=att_probe)
    for row in integrase_attach_audit:
        row["sample"] = args.sample
    audit_rows.extend(integrase_attach_audit)
    rows, candidate_audit = build_candidates(
        args.sample, clusters, systems, contig_lengths,
        args.min_element_bp, args.max_element_bp, args.boundary_bp,
        min_ime_element_bp=args.min_ime_element_bp,
    )
    audit_rows.extend(candidate_audit)

    # --- Phase 3: resolve the real element boundaries -----------------------
    # is_intervals_by_contig was read above, next to the att probe.
    rows, boundary_audit = refine_candidate_boundaries(
        args.sample, rows, args.genome, args.bakta_gff, is_intervals_by_contig,
        args.att_flank_window_bp, args.min_element_bp, args.max_element_bp,
        args.boundary_bp,
    )
    audit_rows.extend(boundary_audit)

    # --- Phase 4b: an ICE cannot live on a plasmid ---------------------------
    # Spec §2.4 is definitional: an ICE INTEGRATES INTO THE CHROMOSOME and encodes
    # conjugation machinery. A plasmid is already its own replicon - it has
    # nothing to integrate into - so a plasmid carrying an integrase, a relaxase
    # and a mating apparatus is a CONJUGATIVE PLASMID that happens to have an
    # integrase, not an ICE.
    #
    # The classifier never knew about replicons; it only appeared to, because
    # plasmid elements happened to lack a clustered integrase. Widening the
    # integrase search window removed that accident and promoted two genuine
    # plasmid elements to "ICE" on these isolates, which is what this pass fixes.
    #
    # The evidence is NOT discarded: the element keeps its anchors and is
    # reported as a conjugative region, with the reason recorded.
    replicon_calls = read_replicon_calls(args.replicons)
    if replicon_calls:
        for row in rows:
            if row.get("element_type") != "ice":
                continue
            call = replicon_calls.get(row.get("contig", ""), "unknown")
            if call == "plasmid":
                row["element_type"] = "conjugative_region"
                row["mge_class"] = "conjugative_region"
                row["mobility"] = "conjugative plasmid region, not an ICE"
                row["mge_id"] = row["mge_id"].replace("|ice-", "|conjugative_region-")
                audit_rows.append(audit_row(
                    args.sample, "kept_flagged", "ice_demoted_on_plasmid_replicon",
                    f"{row['mge_id']}: integrase, relaxase and mating apparatus are "
                    "all present, but this contig is called a PLASMID. An ICE "
                    "integrates into the chromosome (spec §2.4); a plasmid has "
                    "nothing to integrate into, so this is a conjugative plasmid "
                    "carrying an integrase, not an ICE. The anchors are unchanged - "
                    "only the class.",
                    contig=row.get("contig"), start=row.get("start"), end=row.get("end")))
            elif call != "chromosome":
                audit_rows.append(audit_row(
                    args.sample, "kept_flagged", "ice_on_unclassified_replicon",
                    f"{row['mge_id']}: called an ICE, but this contig is not "
                    "classified as chromosome or plasmid, so the definitional "
                    "check (an ICE integrates into a chromosome) could not be "
                    "applied. The call stands; treat it with the caution any "
                    "unknown replicon deserves.",
                    contig=row.get("contig"), start=row.get("start"), end=row.get("end")))

    # --- Phase 6: settle the confidence now the boundaries are known ---------
    rows, confidence_audit = finalise_confidence(
        args.sample, rows, require_trna_boundary=args.require_trna_boundary_for_high)
    audit_rows.extend(confidence_audit)

    return finish(0)


if __name__ == "__main__":
    sys.exit(main())
