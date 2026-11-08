#!/usr/bin/env python3
"""Call a strain's trait profile: which plant-growth, biocontrol, stress-tolerance and
bioremediation traits its genome supports, from sets of genes, with the genes behind
every call and a reason for every decision.

Runs as rule trait_profile (shared/40_annotation.smk), the last per-sample step of the
trait layer, which runs only with traits.run: true. By then each protein Bakta
predicted can carry three kinds of evidence:

  S  structure: the protein's 3Di sequence, predicted from its amino acids by
     ProstT5, searched with Foldseek against EcoFoldDB (Ghaly et al. 2025), a
     curated database of protein structures with known ecological functions:
     nutrient cycling, trace-gas oxidation, plant-microbe interactions, osmotic
     stress tolerance. 3Di is Foldseek's structural alphabet: one letter per
     residue, for how it sits against the residue nearest to it in the folded
     protein. ProstT5, a protein language model, writes these letters straight
     from the amino-acid sequence, with no 3D model built, and Foldseek compares
     the strings much as BLAST compares sequences. EcoFoldDB's catalogue entries
     (842 in version 2.0) sit among about 204,000 background structures that the
     catalogue does not describe;
  Q  the Bakta name: product, gene symbol, KEGG ortholog (KO) and EC number;
  F  the Baktfold name, for a protein Bakta could only call "hypothetical protein"
     and Baktfold named from its structure (only when baktfold.run is on).

Each trait is a row of trait_rules_v1.tsv: one or more alternative sets of required
genes, sometimes with the condition that some of them are co-localised, that is,
lie close together on one contig (--max-gene-gap below). For most traits a set
holds several genes that must all be found, often side by side as in an operon;
for a few, one marker gene is enough, and for some of those genes one structure
hit alone can make the call. Each gene is a row of trait_genes_v1.tsv, the gene
dictionary, saying how it is recognised: by which EcoFoldDB entries, KOs, EC
numbers, Bakta symbols and product-name patterns, which product names veto it, and
whether a structure hit alone is enough (s_alone). This script applies the two
tables to one genome and names no gene or trait of its own: changing what a trait
needs, or how a gene is recognised, means changing a table, where the reason is
written down.

What it reads
-------------
  --bakta-tsv       Bakta's {sample}.tsv from 04.annotation/bakta/{sample}/ (rule
                    annotation). Its cds and sorf rows are the proteins: the same
                    set as Bakta's .faa, which is what the structure search took.
                    An sORF (short open reading frame) is a protein shorter than
                    30 amino acids, too short for Bakta's gene caller; Bakta finds
                    it with a search of its own and keeps it only when it matches
                    a known protein.
  --bakta-faa       Bakta's {sample}.faa, read only for protein lengths, so the
                    proteins too long for the structure search can be listed.
  --hits            the Foldseek table of rule trait_structure_search: the hits of
                    this sample's proteins against the whole EcoFoldDB database,
                    within Foldseek's default limits (e-value up to 10, at most
                    1000 hits per protein); columns query, target, evalue, qcov,
                    tcov, no header. About 3 million rows per genome, so it is
                    streamed, never loaded whole. Plain text or gzip (a name
                    ending in .gz).
  --catalogue       EcoFoldDB_descriptions.txt from the EcoFoldDB folder (rule
                    trait_db or trait_db_local): the catalogue entries, their Gene
                    and their Pathway/Activity. A target missing from it is one of
                    the database's background structures, not a catalogue entry.
  --rules, --genes  the two trait tables, beside this script by default.
  --antismash-json  optional: antiSMASH's {sample}.json from
                    04.annotation/antismash/{sample}/. Its region types are
                    supporting evidence only.
  --baktfold-tsv    optional: Baktfold's {sample}.tsv from
                    04.annotation/baktfold/{sample}/, passed only when
                    baktfold.run is on.
  --max-gene-gap    traits.max_gene_gap: how far apart, in genes along a contig,
                    each gene of a co-localised group may lie from the one before
                    it (1 = next-door genes only).
  --max-protein-length  the length limit of the structure search (proteins longer
                    than this were not searched).

What it writes
--------------
  --out-genes    {sample}_trait_genes.tsv: one row per protein and dictionary gene
                 with any evidence, counted or not, and why.
  --out-profile  {sample}_trait_profile.tsv: one row per trait, in table order:
                 called, partial or not called, with the genes behind it.
                 Rule traits_summary reads every sample's profile.
  --out-audit    {sample}_trait_profile_audit.tsv: one row per trait with a fixed
                 reason token and the detail behind the decision.
  --out-provenance  {sample}_trait_provenance.json (optional; the workflow always
                 asks for it): the SHA-256 fingerprints of the two trait tables,
                 the EcoFoldDB catalogue and the three tables above, with the
                 calling parameters. Rules traits_summary and traits_workbook
                 check them before quoting the catalogue beside these calls; see
                 write_provenance().
Column meanings and every reason token are listed beside GENE_TABLE_COLUMNS,
PROFILE_COLUMNS and AUDIT_COLUMNS below, and in docs/analysis/traits.md.

How a protein counts for a gene
-------------------------------
  1. Structure: the protein's single best hit in the whole Foldseek table (lowest
     e-value; the first listed wins a tie) must be a catalogue entry with e-value
     below 1e-10 and query and target coverage above 0.8, EcoFoldDB's published
     filter. The e-value is the number of hits this good expected by chance in a
     database this size, so lower is better; query and target coverage are the
     shares of our protein and of the database structure that the alignment
     spans. A kept hit gives S evidence to every dictionary gene whose ecofolddb
     list holds the entry's key (Gene@Pathway/Activity).
  2. Names: a product_exclude pattern matching the product vetoes the gene; else a
     product_include pattern, a KO, an EC number or the exact Bakta symbol matches
     (a symbol marked symbol_needs_product counts only with the product or a KO).
     A product is read without the prefix Bakta writes before the name of a
     protein it marks: "(pseudo) " for a pseudogene, "(5' truncated) ",
     "(3' truncated) " or "(partial) " for a gene that runs off an end of its
     sequence (point 5). Baktfold's products carry the same prefixes and are read
     the same way.
  3. A veto, from the Bakta product or the Baktfold product, cancels every kind of
     evidence for that gene on that protein, structure included, even for a gene
     whose structure hit would count alone. Why structure too: related enzymes
     often share one fold, so a structure hit can land on the wrong member of a
     family, while the name usually says which member the protein is.
     product_exclude lists those relatives, such as MiaA, which prenylates tRNA,
     for the cytokinin enzyme ipt.
  4. The gene counts on the protein when it is not vetoed and has Q or F, or has S
     and its s_alone column says yes, unless the protein is a pseudogene (point 5)
     or point 6 gives the protein to another gene. One protein may count for
     several genes, such as a fused enzyme whose name says so.
  5. Bakta's marks. A gene marked as running off the start or end of its sequence
     counts as usual, with a note in the detail. Why: BacFlux runs Bakta without
     --partial, so these marks appear only on circular replicons, which have no
     ends. A pseudogene never counts: aligning it to its reference protein, Bakta
     found the gene broken by a frameshift, an internal stop codon, or a lost
     start codon that moves the start inside the gene (Bakta 1.12,
     features/cds.py), so it may not make a working protein. Where its evidence
     would count, the gene table lists it as not counted with reason pseudogene,
     and the trait calls take the protein as a broken copy of the gene, unless
     point 6 gives it to another gene; a set that is complete except for
     pseudogenes is reported as partial (below). From here on, a gene is found
     as a pseudogene when the trait calls take a pseudogene as it, and found
     only as a pseudogene when, in addition, no protein counts for it.
  6. One protein, one gene per rule. A protein counts for two genes that one rule
     names (in its required sets or its supporting list) only when one source
     names both: one Bakta name (a fusion such as PqqC/D, whose product matches
     both genes, or a product that matches one gene and a KO or EC number that
     matches the other), one Baktfold name, or one structure hit whose catalogue
     entry both genes list. When the two genes come from different sources, so
     that their evidence letters have none in common, the protein counts once:
     for the structure's gene when a structure hit counts alone for that gene
     (s_alone yes), even against a Bakta and a Baktfold name together;
     otherwise the names decide, and the gene of Bakta's name ranks before that
     of Baktfold's, even when the structure hit points to the gene of
     Baktfold's name. The other gene keeps its row, not counted, with reason
     counted_for_other_gene; this holds wherever that gene is used, so such a
     protein never fills two places of a set. Genes that share no rule are not
     compared, and rows that do not count anyway (vetoed, or structure alone for
     a gene with s_alone no) take no part. A pseudogene is compared in the same
     way, on its rows whose evidence would count (reason pseudogene): the trait
     calls take the broken protein as the gene that ranks first, and the other
     gene's row keeps reason pseudogene, its detail naming the gene the protein
     is taken as. So such a broken protein never fills two places of a set
     either. Why: a protein is one gene product. A structure hit to the alpha
     chain of protocatechuate 3,4-dioxygenase (PcaG) on a protein Bakta names as
     the beta chain (PcaH) must not fill both places of the pcaG+pcaH set. Here
     the protein counts for pcaG, because a structure hit counts alone for pcaG
     (s_alone yes). And why the structure hit decides only for such a gene: a
     gene has s_alone no because related enzymes share its fold, so a structure
     hit can land on the wrong member of the family (point 3). Baktfold's name
     rests on the same predicted structure as the hit (ProstT5's 3Di string,
     searched with Foldseek), while Bakta's name comes from sequence similarity,
     which usually tells such relatives apart.

How a trait is decided
----------------------
  called      the first required set (in table order) whose genes are all counted
              and whose co-localised groups each lie in one chain on one contig:
              each gene at most --max-gene-gap genes from the one before it;
  partial     otherwise, the first set (in table order) that is one of these and
              is not withheld:
                - two or more genes, all counted but not co-localised;
                - every gene counted or a pseudogene, at least one a pseudogene;
                  this holds for a set of one gene too (its gene is there, but
                  broken), which is otherwise never partial;
                - two or more genes, at least half of them counted (rounded up).
              A partial of the first or the last kind is withheld when the rule
              names key_genes, the genes that set the pathway apart from
              look-alikes, and none of them is counted, or found as a pseudogene,
              anywhere. One of the middle kind never is: every gene of that set is
              in the genome, and the partial is there to show the broken one;
  not called  otherwise, including a trait whose partials were all withheld.
Supporting genes and antiSMASH regions are reported beside a call and never count
toward it. A supporting gene can still change a call through point 6: a protein
that counts for it does not count for a required gene of the same rule that a
different source names on that protein.

In the audit detail and the profile's missing column, a gene that is not counted is
written as gene:locus_tag(pseudogene) when it is found only as a pseudogene
(point 5), and as its bare name otherwise: when no protein matched it, or none that
counts for it or is taken as it (a vetoed protein, say, or one that point 6 gives
to another gene). The audit detail of a partial, or of a trait that is not called,
names every required gene of the trait found only as a pseudogene, those outside
the set it reports included, so that, in a trait that is partial or not called, no
such gene shows in the gene table alone. A pseudogene that point 6 gives to another
gene is never written under the gene that lost it: that gene's row in the gene
table keeps reason pseudogene and names the gene the protein is taken as.

What a call does not mean: the genes are present, not that they are expressed or
that the trait works in this strain. The caveat column of every profile row says
what that call leaves open, and the assay column of trait_rules_v1.tsv names the
laboratory test that would confirm it.

Standard library only (Python 3.12 or later). Rule trait_profile runs it in
workflow/envs/python.yaml.
"""

import argparse
import collections
import decimal
import gzip
import hashlib
import json
import math
import os
import re
import sys


# ── The structure filter: EcoFoldDB's published one ──────────────────────────
# A protein's structure evidence is its single best Foldseek hit across the whole
# EcoFoldDB database. It counts only when that hit is a catalogue entry with an
# e-value below E_VALUE_MAX and both query and target coverage above COVERAGE_MIN
# (strict inequalities). These are the filters the EcoFoldDB authors published with
# the database (Ghaly et al. 2025, Environ Microbiol 27:e70178), so a kept hit here
# means what it means in their paper. Taking the best hit across the WHOLE database,
# background structures included, is what makes the filter strict: a protein whose
# closest structure is not a catalogue entry gets no structure evidence at all.
E_VALUE_MAX = 1e-10
COVERAGE_MIN = 0.8

# ── Borderline: a kept hit close to a cut-off or to a rival ──────────────────
# A small change in a predicted structure can move a hit across a threshold, so a
# kept hit is flagged borderline, without changing any decision, when:
#   - its e-value is within BORDERLINE_FACTOR of the cut-off (e >= 1e-11);
#   - its query or target coverage is at or below BORDERLINE_COVERAGE;
#   - the best hit to a DIFFERENT catalogue gene, or to a target outside the
#     catalogue, has an e-value within BORDERLINE_FACTOR of the best hit's.
# The flag tells the reader which structure calls deserve a second look. The values
# are a starting point, not a calibrated cut-off.
BORDERLINE_FACTOR = 10
BORDERLINE_COVERAGE = 0.85
# E_VALUE_MAX / BORDERLINE_FACTOR, written out: the division in floating point
# gives 1.0000000000000001e-11, and a hit Foldseek prints as 1.000E-11 would then
# miss the flag it should get.
BORDERLINE_E_VALUE = 1e-11

# ── Defaults for a hand-run ───────────────────────────────────────────────────
# Rule trait_profile passes both numbers: --max-gene-gap from traits.max_gene_gap,
# and --max-protein-length from TRAITS_MAX_PROTEIN_LENGTH, the same limit rule
# trait_structure_query applies before ProstT5 (EcoFoldDB's published method).
# The defaults here are what a hand-run of the script gets.
DEFAULT_MAX_GENE_GAP = 10
DEFAULT_MAX_PROTEIN_LENGTH = 4000

# ── The two trait tables shipped beside this script ──────────────────────────
# Found through this file's own location, so a hand-run works from any folder.
# Rule trait_profile passes them explicitly (TRAIT_RULES_TSV, TRAIT_GENES_TSV).
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_RULES_TSV = os.path.join(_HERE, "trait_rules_v1.tsv")
DEFAULT_GENES_TSV = os.path.join(_HERE, "trait_genes_v1.tsv")

# The columns each table must have, found by name in its header line. What each
# column holds: docs/analysis/traits.md; how this script reads them:
# load_trait_tables.
RULE_COLUMNS = ["trait_id", "trait", "application", "required", "supporting",
                "colocalised", "key_genes", "caveat", "assay", "reference"]
GENE_COLUMNS = ["gene", "protein", "ecofolddb", "ko", "ec", "bakta_gene",
                "symbol_needs_product", "product_include", "product_exclude",
                "s_alone", "notes"]

# A supporting token with this prefix names an antiSMASH region type, not a gene,
# for example antiSMASH:NRP-metallophore.
ANTISMASH_PREFIX = "antiSMASH:"

# ── What the other inputs give us ────────────────────────────────────────────
# The columns read from Bakta's annotation table, found by name in the table's own
# header line ("#Sequence Id<TAB>Type<TAB>Start..."). Baktfold writes its
# {sample}.tsv in the same layout, so one reader serves both.
BAKTA_COLUMNS_NEEDED = ["Sequence Id", "Type", "Start", "Stop", "Strand",
                        "Locus Tag", "Gene", "Product", "DbXrefs"]

# The feature types that are proteins in Bakta's .faa, and so were given to the
# structure search: ordinary CDS and short open reading frames (sORFs, proteins
# shorter than 30 amino acids; module docstring, --bakta-tsv).
PROTEIN_TYPES = {"cds", "sorf"}

# The prefixes Bakta writes before the product of a protein it marks, in the
# Product column of {sample}.tsv (Bakta 1.12, function write_features), and the
# mark each one stands for. A protein carries at most one, "(pseudo) " first:
#   "(pseudo) "         a pseudogene: aligned to its reference protein, the gene is
#                       broken by a frameshift, an internal stop codon, or a lost
#                       start codon that moves the start inside the gene;
#   "(5' truncated) "   the gene runs off an end of its sequence at its 5' end,
#                       where translation starts;
#   "(3' truncated) "   the same at its 3' end, where translation stops;
#   "(partial) "        the gene runs off both ends of its sequence.
# BacFlux runs Bakta without --partial, so on a linear contig Bakta's gene caller
# works with closed ends: it predicts only genes whose start and stop codons both
# lie inside the contig, and none carries the last three marks (Bakta 1.12,
# features/cds.py). They appear only on circular replicons (nanopore and hybrid
# modes), most often on the replication gene at position 1, such as dnaA: the
# long-read modes turn each circular contig with dnaapler so that it starts there.
# A circular replicon has no ends, so such a gene is not cut in the genome, and it
# counts as usual.
# Bakta's .faa keeps all these proteins, so they were searched by structure too.
# read_annotation_table splits the prefix off, and every name test reads the
# product without it. Why: a pattern anchored with "^" never matches a product with
# a prefix, while an unanchored pattern or a KO does, so with the prefix left on, a
# pseudogene would count through some genes and not through others. Baktfold
# writes its .tsv in Bakta's layout, with the same prefixes, and is read the same
# way. What each mark does to the evidence is decided in collect_evidence.
BAKTA_PRODUCT_PREFIXES = {
    "(pseudo) ": "pseudogene",
    "(5' truncated) ": "5' truncated",
    "(3' truncated) ": "3' truncated",
    "(partial) ": "partial",
}

# The few words each mark adds to the detail of the protein's gene-table rows. A
# note names the end of the gene, not of the sequence: Bakta writes "(5' truncated)"
# for a forward-strand gene that runs off the start of its sequence and for a
# reverse-strand gene that runs off its end (Bakta 1.12, features/cds.py).
BAKTA_MARK_NOTES = {
    "pseudogene": "Bakta marks the gene as a pseudogene (disrupted reading frame)",
    "5' truncated": "Bakta: runs off its sequence at the gene's 5' end (5' truncated)",
    "3' truncated": "Bakta: runs off its sequence at the gene's 3' end (3' truncated)",
    "partial": "Bakta: runs off both ends of its sequence (partial)",
}

# Bakta's product for a protein it could not name. Only a protein with exactly this
# product, read without Bakta's prefix, can carry a Baktfold name (F evidence). As
# BacFlux runs it, Baktfold renames only the CDS that Bakta flagged hypothetical,
# pseudogenes included ("(pseudo) hypothetical protein"). Such a pseudogene takes
# its Baktfold name like any other protein, so the gene table shows what Baktfold
# saw in it, and collect_evidence then keeps it from counting, as it does for
# every pseudogene.
HYPOTHETICAL_PRODUCT = "hypothetical protein"

# The catalogue columns read from EcoFoldDB_descriptions.txt.
CATALOGUE_COLUMNS_NEEDED = ["ID", "Gene", "Pathway/Activity"]

# The columns of the Foldseek table, in this order. Rule trait_structure_search asks
# Foldseek's convertalis step, which writes the search results out as a table, for
# exactly these (--format-output query,target,evalue,qcov,tcov), so the two must be
# changed together.
FOLDSEEK_COLUMNS = ["query", "target", "evalue", "qcov", "tcov"]


# ── What we write: three tables and their vocabulary ─────────────────────────
# Every table starts with the sample, so tables of many genomes can be stacked.
# "NA" marks a cell with nothing to show. Lists inside a cell are joined by ", "
# unless said otherwise.
#
# {sample}_trait_genes.tsv, one row per protein and dictionary gene with any
# evidence, plus the rows about a protein rather than a gene (gene NA). Sorted by
# contig (in the order Bakta lists them), cds_index, gene, with a protein's own
# rows (gene NA) after its gene rows.
#   sample            the --sample argument
#   locus_tag         Bakta's locus tag of the protein
#   contig, start, stop, strand
#                     where Bakta placed the gene: 1-based and inclusive; a start
#                     larger than the stop means the gene runs across the end of a
#                     circular contig
#   cds_index         the protein's rank by start among this contig's proteins (cds
#                     and sORF), counting from 0; co-localisation is measured in
#                     these steps
#   gene              the dictionary gene (trait_genes_v1.tsv), or NA
#   traits            the trait_ids whose required or supporting genes name it
#   counted           yes when the protein counts for the gene in the trait calls
#   evidence          what was seen, joined by "+" in this order: S structure, Q
#                     Bakta name, F Baktfold name. On a vetoed row: what the veto
#                     cancelled
#   bakta_gene, bakta_product, ko, ec
#                     Bakta's symbol, product, KEGG orthologs and EC numbers. The
#                     product is as Bakta wrote it, with the prefix of a protein
#                     Bakta marks ("(pseudo) ", "(5' truncated) ",
#                     "(3' truncated) ", "(partial) "); every name test reads it
#                     without the prefix
#   baktfold_product  Baktfold's name, as Baktfold wrote it, for a protein Bakta
#                     called "hypothetical protein" and Baktfold named otherwise
#                     (prefixes aside); NA for every other protein
#   ecofolddb_target  the protein's best Foldseek hit in EcoFoldDB, on every row of
#                     a protein that has one, whether or not it counted
#   ecofolddb_gene, ecofolddb_pathway
#                     the catalogue's Gene and Pathway/Activity of that target; NA
#                     when the target is outside the catalogue
#   evalue, qcov, tcov
#                     that hit's e-value, exactly as Foldseek printed it, and its
#                     query and target coverage
#   borderline        on a row whose evidence includes a kept structure hit (S):
#                     the borderline reasons joined by "; ", or "no"; NA on other
#                     rows
#   decision          counted | not counted | not searched
#   reason            a fixed token, listed below
#   detail            what matched, or why it did not count, in plain words. When
#                     point 6 of the module docstring gives a protein to one gene
#                     rather than another (two genes that share a rule, and no one
#                     source names both), the row of the gene it is given to names
#                     the other gene, and the other gene's row
#                     (counted_for_other_gene, or pseudogene on a broken protein)
#                     names the gene it is given to and why. On a gene row of a
#                     protein Bakta marks as a pseudogene or as running off an end
#                     of its sequence, it ends with a note saying so
GENE_TABLE_COLUMNS = [
    "sample", "locus_tag", "contig", "start", "stop", "strand", "cds_index",
    "gene", "traits", "counted", "evidence",
    "bakta_gene", "bakta_product", "ko", "ec", "baktfold_product",
    "ecofolddb_target", "ecofolddb_gene", "ecofolddb_pathway",
    "evalue", "qcov", "tcov", "borderline",
    "decision", "reason", "detail",
]

# {sample}_trait_profile.tsv, one row per trait, in the order of trait_rules_v1.tsv.
#   sample
#   trait_id, trait, application
#                     copied from trait_rules_v1.tsv
#   decision          called | partial | not called
#   genes             called or partial: the counted genes of that set, each as
#                     gene:locus_tag(evidence), several proteins of one gene joined
#                     by ","; genes joined by "; ". NA when not called, and on a
#                     partial none of whose genes is counted
#   missing           partial: the genes of that set not counted, joined by ", ".
#                     A gene found only as a pseudogene (module docstring,
#                     point 5) is written gene:locus_tag(pseudogene), several
#                     proteins joined by ","; any other gene not counted, its
#                     bare name. NA otherwise
#   supporting        the counted supporting genes, in the same form, and each
#                     antiSMASH region type found over a protein counted for a
#                     required gene, as antiSMASH:<type>(contig:start-end), 1-based
#                     and inclusive like Bakta's (a start larger than the end means
#                     the region runs across the end of a circular contig), or
#                     antiSMASH:<type>(not checked) when no antiSMASH output was
#                     given. Reported whatever the decision, and never counted
#                     toward it (but see module docstring, point 6)
#   borderline        called or partial: yes when a protein counted for the call
#                     carries a borderline structure hit for its gene, else no.
#                     NA when not called
#   structure_only    called or partial: the genes whose every counted protein
#                     rests on structure alone (evidence S), with no name
#                     agreeing. NA when there are none
#   caveat            copied from trait_rules_v1.tsv: what the call does not show
PROFILE_COLUMNS = [
    "sample", "trait_id", "trait", "application", "decision", "genes",
    "missing", "supporting", "borderline", "structure_only", "caveat",
]

# {sample}_trait_profile_audit.tsv, one row per trait: the decision, a fixed reason
# token and the detail behind it (the set, the genes found and missing).
AUDIT_COLUMNS = ["sample", "trait_id", "decision", "reason", "detail"]

# EVERY decision / reason PAIR THIS SCRIPT CAN WRITE.
#
# {sample}_trait_genes.tsv
#   decision=counted
#     name_and_structure_agree   a Bakta or Baktfold name and the kept structure
#                                hit both point to the gene
#     name_only                  a name points to the gene and the structure does
#                                not; the detail says what the best hit was
#     structure_counts_alone     only the kept structure hit points to the gene,
#                                and this gene's s_alone column says yes
#   decision=not counted
#     pseudogene                 the evidence would count the gene (one of the three
#                                reasons above), but Bakta marks the protein as a
#                                pseudogene ("(pseudo) " before its product); the
#                                detail says what matched. The trait calls take the
#                                protein as a broken copy of the gene on the same
#                                terms as a working protein counts (module
#                                docstring, points 5 and 6): when point 6 gives it
#                                instead to another gene of one of this gene's
#                                rules, named by a different source that ranks
#                                first, this row keeps its reason, and its detail
#                                names that gene
#     counted_for_other_gene     the evidence would count the gene, but the protein
#                                counts instead for another gene of one of its
#                                rules, named by a different source that ranks
#                                first: the structure hit, when it counts alone
#                                for its gene (s_alone yes), then Bakta's name,
#                                then Baktfold's name (module docstring, point 6:
#                                one protein, one gene per rule); the detail
#                                names both genes and the evidence for each
#     structure_only_needs_name  only the kept structure hit points to the gene,
#                                and this gene counts from structure only when a
#                                name agrees (s_alone no)
#     vetoed_by_name             the gene's product_exclude pattern matches the
#                                protein's Bakta or Baktfold product, so none of
#                                its evidence counts, structure included
#     catalogue_gene_unused      (gene NA) a kept structure hit to a catalogue entry
#                                that no dictionary gene uses
#   decision=not searched
#     not_searched_too_long      (gene NA) the protein is longer than the structure
#                                search takes; its names were still checked
# A gene Bakta marks as running off an end of its sequence has no reason of its
# own: it counts as usual, and the detail ends with a note naming the mark.
#
# {sample}_trait_profile_audit.tsv
#   decision=called
#     set_complete               every gene of a required set found, and
#                                co-localised wherever the trait asks
#   decision=partial
#     set_complete_not_colocalised  every gene of a set of two or more found, but
#                                a group that must lie together does not
#     set_complete_but_pseudogene  every gene of a set found, counted or as a
#                                pseudogene, and at least one only as a
#                                pseudogene; a set of one gene can be partial this
#                                way. The detail names each broken gene with its
#                                locus tag
#     half_of_set_found          at least half (rounded up) of a set of two or more
#                                genes found
#   decision=not called
#     partial_withheld           a set_complete_not_colocalised or half_of_set_found
#                                partial would apply, but none of the trait's
#                                key_genes is found anywhere, counted or as a
#                                pseudogene, and no set qualifies as
#                                set_complete_but_pseudogene, which key_genes never
#                                withholds; the detail reads "no key gene (desA)
#                                found (...)", with the partial it withheld
#     too_few_genes              some genes found, fewer than a partial needs; the
#                                detail names the closest set
#     no_gene_found              no gene of any required set counted; the detail
#                                names any required gene found only as a
#                                pseudogene
# Here "found" means counted, except where a pseudogene is named ("found as a
# pseudogene" and "found only as a pseudogene": module docstring, point 5). A
# detail that lists genes not counted writes a gene found only as a pseudogene as
# gene:locus_tag(pseudogene), and any other gene not counted as its bare name,
# including a gene whose pseudogenes are all given to another gene (point 6). The
# detail of a partial, of partial_withheld and of too_few_genes reports one set; it
# ends by naming the trait's required genes outside that set found only as a
# pseudogene, as "; outside this set, found only as a pseudogene:
# gene:locus_tag(pseudogene), ...". no_gene_found names them in its own words.


# ── Reading Bakta's (and Baktfold's) annotation table ────────────────────────
# Bakta's {sample}.tsv lists every feature it annotated; its cds and sorf rows are
# the proteins the trait profile works on. Baktfold writes its {sample}.tsv in the
# same layout, so one reader serves both. This section gives collect_evidence the
# proteins with their names, Bakta's marks and their places on the contig
# (read_annotation_table, add_positions), the Baktfold names (read_baktfold_names)
# and the protein lengths (read_protein_lengths).

def _int_or_none(value):
    """Parse an integer field, returning None instead of raising when the field
    is missing or not a number. read_annotation_table uses it for Bakta's Start
    and Stop, and stops the run on a None there.

    Copied from scripts/80_mobilome/isescan_to_table.py rather than imported:
    Python looks for modules in the running script's own folder, and each folder
    of workflow/scripts belongs to one rule file (its README.md).
    """
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def parse_dbxrefs(text):
    """Pull the KEGG orthologs and EC numbers out of Bakta's DbXrefs cell.

    Takes in: the DbXrefs cell, e.g. "EC:3.5.1.5, KEGG:K01428, SO:0001217".
    Returns:  (list of KOs, list of EC numbers), in the order written, e.g.
              (["K01428"], ["3.5.1.5"]). Every other cross-reference is ignored.
    """
    kos = []
    ecs = []
    for item in text.split(","):
        item = item.strip()
        if item.startswith("KEGG:"):
            kos.append(item[len("KEGG:"):])
        elif item.startswith("EC:"):
            ecs.append(item[len("EC:"):])
    return kos, ecs


def split_bakta_prefix(product):
    """Split the prefix Bakta writes before the product of a protein it marks.

    Takes in: one Product cell of Bakta's (or Baktfold's) .tsv, for example
              "(pseudo) Alkaline phosphatase".
    Returns:  (the product without its prefix, the mark), for example
              ("Alkaline phosphatase", "pseudogene"). A product with none of the
              prefixes of BAKTA_PRODUCT_PREFIXES comes back whole, with mark None.
    """
    for prefix, mark in BAKTA_PRODUCT_PREFIXES.items():
        if product.startswith(prefix):
            return product[len(prefix):], mark
    return product, None


def read_annotation_table(path):
    """Read the proteins of a Bakta (or Baktfold) {sample}.tsv, in file order.

    Takes in: the table. It opens with a few comment lines starting with '#'; the
              last of them is the column header (#Sequence Id, Type, Start, Stop,
              Strand, Locus Tag, Gene, Product, DbXrefs).
    Does:     keeps the cds and sorf rows, the proteins of Bakta's .faa, and splits
              off the prefix of a marked protein's product (split_bakta_prefix).
    Returns:  a list of dicts {contig, start, stop, strand, locus_tag, gene,
              product, product_written, mark, kos, ecs}. product is the name
              without Bakta's prefix, the one every name test reads;
              product_written is the cell as written, for the output tables;
              mark is what the prefix said ("pseudogene", "5' truncated",
              "3' truncated", "partial"), or None.

    The header search and column lookup are copied from read_bakta_genes in
    scripts/80_mobilome/isescan_to_table.py. One thing differs on purpose: a
    protein row that is too short or has no numeric coordinates stops the run here
    instead of being skipped, because a protein dropped in silence could turn a
    called trait into a partial one. (A row too short even to show its Type says
    nothing about a protein, and is skipped.)
    """
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n").rstrip("\r") for line in handle]

    # The column header is the comment line that starts with '#Sequence Id'.
    header = None
    for line in lines:
        if line.startswith("#Sequence Id"):
            header = line[1:].split("\t")
    if header is None:
        raise ValueError(
            f"{path}: no '#Sequence Id' header line, so this does not look like "
            "Bakta's .tsv annotation table."
        )
    missing = [name for name in BAKTA_COLUMNS_NEEDED if name not in header]
    if missing:
        raise ValueError(
            f"{path}: annotation table is missing column(s) {', '.join(missing)}. "
            f"Header seen: {header}"
        )
    column = {name: header.index(name) for name in BAKTA_COLUMNS_NEEDED}
    last_needed = max(column.values())

    proteins = []
    for number, line in enumerate(lines, start=1):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) <= column["Type"]:
            continue                      # too short to say what it is
        if fields[column["Type"]].strip() not in PROTEIN_TYPES:
            continue                      # tRNA, rRNA, ncRNA, oriC, ...
        if len(fields) <= last_needed:
            raise ValueError(f"{path}: line {number} is a protein row with only "
                             f"{len(fields)} column(s).")
        start = _int_or_none(fields[column["Start"]])
        stop = _int_or_none(fields[column["Stop"]])
        if start is None or stop is None:
            raise ValueError(f"{path}: line {number} is a protein row without "
                             "numeric Start and Stop.")
        kos, ecs = parse_dbxrefs(fields[column["DbXrefs"]])
        product_written = fields[column["Product"]].strip()
        product, mark = split_bakta_prefix(product_written)
        proteins.append({
            "contig": fields[column["Sequence Id"]].strip(),
            "start": start,
            "stop": stop,
            "strand": fields[column["Strand"]].strip(),
            "locus_tag": fields[column["Locus Tag"]].strip(),
            "gene": fields[column["Gene"]].strip(),
            "product": product,
            "product_written": product_written,
            "mark": mark,
            "kos": kos,
            "ecs": ecs,
        })
    return proteins


def add_positions(proteins):
    """Give every protein its contig's order and its CDS index, in place.

    Takes in: the proteins from read_annotation_table, in Bakta's file order.
    Does:     ranks the proteins of each contig by start, counting from 0. Sorting
              is stable, so two proteins with the same start keep Bakta's order.
              Contigs are treated as linear: the first and last gene of a circular
              chromosome are not neighbours here. Why: Bakta's table does not say
              which contigs are circular, and in a draft assembly most contigs are
              fragments anyway. Treating a contig as linear can only break a
              chain, never make one: a cluster that straddles the origin of a
              closed chromosome or plasmid reads as not co-localised, so its trait
              comes out partial (set_complete_not_colocalised), never wrongly
              called.
    Adds:     "cds_index" (the rank) and "contig_order" (the contig's position in
              Bakta's table, used to sort the gene table the way Bakta does).
    """
    contig_order = {}
    by_contig = collections.defaultdict(list)
    for protein in proteins:
        if protein["contig"] not in contig_order:
            contig_order[protein["contig"]] = len(contig_order)
        by_contig[protein["contig"]].append(protein)
    for contig, members in by_contig.items():
        members.sort(key=lambda protein: protein["start"])
        for index, protein in enumerate(members):
            protein["cds_index"] = index
            protein["contig_order"] = contig_order[contig]


def read_baktfold_names(path, proteins):
    """Find the proteins that carry a Baktfold name (F evidence).

    Takes in: Baktfold's {sample}.tsv (rule structure_annotation, only when
              baktfold.run is on) and Bakta's proteins.
    Does:     a protein has a Baktfold name when Bakta's product is exactly
              "hypothetical protein" and Baktfold gave it another product. Both
              products are compared without Bakta's prefix (read_annotation_table
              splits it off both tables), so a pseudogene such as
              "(pseudo) hypothetical protein" can take a Baktfold name too. Its
              Baktfold symbol, product, KOs and EC numbers are kept. Whether a
              pseudogene counts is decided later, in collect_evidence, from the
              mark in Bakta's table.
    Returns:  {locus_tag: {gene, product, product_written, mark, kos, ecs}} for
              those proteins only, as read_annotation_table reads Baktfold's rows.
    """
    baktfold_by_locus = {row["locus_tag"]: row for row in read_annotation_table(path)}
    names = {}
    for protein in proteins:
        # Only a protein Bakta could not name can carry a Baktfold name.
        if protein["product"] != HYPOTHETICAL_PRODUCT:
            continue
        renamed = baktfold_by_locus.get(protein["locus_tag"])
        if renamed is None:
            continue
        # Baktfold found no name either: no F evidence.
        if not renamed["product"] or renamed["product"] == HYPOTHETICAL_PRODUCT:
            continue
        names[protein["locus_tag"]] = renamed
    return names


def read_protein_lengths(path):
    """Read Bakta's .faa and return {protein id: length in amino acids}.

    The id is the first word of the header line, Bakta's locus tag. Sequences may
    be wrapped over several lines; the length is the sum of the line lengths,
    which is what the length filter of rule trait_structure_query counts too.
    The lengths tell which proteins were too long for the structure search:
    collect_evidence lists them as not searched, and main counts the searchable
    ones before warning about an empty hit table.
    """
    lengths = {}
    current = None
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                words = line[1:].split()
                current = words[0] if words else ""
                lengths[current] = 0
            elif current is not None:
                lengths[current] += len(line)
    return lengths


# ── Reading EcoFoldDB: the catalogue and the Foldseek hits ───────────────────
# The catalogue names the gene and the pathway of each catalogue entry; the
# Foldseek table lists the hits of every searched protein. read_best_hits keeps
# each protein's best hit and its best rival, and judge_best_hits applies
# EcoFoldDB's filter and the borderline flags. The result is the S evidence that
# collect_evidence reads.

def read_catalogue(path):
    """Read EcoFoldDB_descriptions.txt into {target ID: {gene, pathway, key}}.

    Takes in: the descriptions table of the EcoFoldDB folder, one row per
              catalogue entry (ID, Gene, Pathway/Activity and other columns).
    Returns:  the catalogue. A Foldseek target missing from it is a background
              structure. The key, Gene + "@" + Pathway/Activity, is how
              trait_genes_v1.tsv names catalogue entries: one key covers every
              entry with that Gene and Pathway/Activity (sometimes several
              structures), and the same Gene in another pathway has another key.
    Used by read_best_hits and judge_best_hits (what a hit points to) and by
    load_trait_tables (every ecofolddb key must be in it).
    """
    catalogue = {}
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").rstrip("\r").split("\t")
        missing = [name for name in CATALOGUE_COLUMNS_NEEDED if name not in header]
        if missing:
            raise ValueError(
                f"{path}: missing column(s) {', '.join(missing)}, so this does not "
                f"look like EcoFoldDB_descriptions.txt. Header seen: {header}"
            )
        id_index = header.index("ID")
        gene_index = header.index("Gene")
        pathway_index = header.index("Pathway/Activity")
        last_needed = max(id_index, gene_index, pathway_index)
        for line in handle:
            fields = line.rstrip("\n").rstrip("\r").split("\t")
            if len(fields) <= last_needed:
                continue                  # a blank or broken line names no entry
            gene = fields[gene_index].strip()
            pathway = fields[pathway_index].strip()
            catalogue[fields[id_index].strip()] = {
                "gene": gene,
                "pathway": pathway,
                "key": gene + "@" + pathway,
            }
    return catalogue


def open_text(path):
    """Open the Foldseek hit table (--hits) for reading as text: plain, or
    gzip-compressed when the name ends in .gz. read_best_hits streams the table
    through it, line by line. Rule trait_profile passes the plain table of rule
    trait_structure_search; a gzip-compressed copy, for a hand-run, reads the same.

    Copied from scripts/80_mobilome/check_gtdb_organism_table.py, which streams
    the GTDB metadata the same way.
    """
    opener = gzip.open if path.endswith(".gz") else open
    return opener(path, "rt", newline="", encoding="utf-8", errors="replace")


# One Foldseek hit, kept only while it is a protein's best hit or its best rival.
# evalue is the number hits are compared by; label is the catalogue Gene of the
# target, or None for a target outside the catalogue. The texts are kept as
# Foldseek printed them, for the output and for the near-tie test.
Hit = collections.namedtuple(
    "Hit", ["evalue", "label", "target", "evalue_text", "qcov_text", "tcov_text"])


def read_best_hits(path, catalogue):
    """Stream the Foldseek table and keep, per protein, its best hit and best rival.

    Takes in: the hit table of rule trait_structure_search (FOLDSEEK_COLUMNS, no
              header) and the catalogue.
    Does:     the best hit is the one with the lowest e-value; with a strict "<",
              the first listed wins a tie. The rival is the best hit whose label
              differs from the best hit's: another catalogue Gene, or a target
              outside the catalogue (all such targets share one label, None).
              Both are updated hit by hit, so memory grows with the number of
              proteins, not with the about 3 million rows. When a new best hit
              arrives with a different label, the old best becomes the rival: it
              was better than every earlier hit, so it is the best of the other
              labels. When the new best has the same label, the rival stands.
    Returns:  ({query: best Hit}, {query: rival Hit}); a protein with no rival is
              absent from the second dict.

    A line that cannot be read stops the run: a damaged hit table is a broken
    input, and carrying on would drop structure evidence without a trace.
    """
    best = {}
    rival = {}
    with open_text(path) as handle:
        for number, line in enumerate(handle, start=1):
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < len(FOLDSEEK_COLUMNS):
                raise ValueError(
                    f"{path}: line {number} has {len(fields)} column(s); expected "
                    f"{', '.join(FOLDSEEK_COLUMNS)} (foldseek convertalis)."
                )
            query, target, evalue_text, qcov_text, tcov_text = fields[:5]
            try:
                evalue = float(evalue_text)
            except ValueError:
                raise ValueError(
                    f"{path}: line {number}: e-value '{evalue_text}' is not a number."
                ) from None
            # float() also accepts 'nan' and 'inf'. Foldseek never prints them, and a
            # 'nan' would slip through every comparison below, so it stops the run too.
            if not math.isfinite(evalue):
                raise ValueError(
                    f"{path}: line {number}: e-value '{evalue_text}' is not a finite "
                    "number."
                )
            entry = catalogue.get(target)
            label = entry["gene"] if entry else None

            current = best.get(query)
            if current is None:
                # The protein's first hit is its best so far.
                best[query] = Hit(evalue, label, target, evalue_text, qcov_text, tcov_text)
            elif evalue < current.evalue:
                # A better hit: the old best becomes the rival when its label
                # differs, and the new hit is the best.
                if label != current.label:
                    rival[query] = current
                best[query] = Hit(evalue, label, target, evalue_text, qcov_text, tcov_text)
            elif label != current.label:
                # Not better than the best, but of another label: it becomes the
                # rival when it beats the rival so far.
                previous_rival = rival.get(query)
                if previous_rival is None or evalue < previous_rival.evalue:
                    rival[query] = Hit(evalue, label, target, evalue_text, qcov_text,
                                       tcov_text)
    return best, rival


def _coverage(text, query, path):
    """Parse a coverage value of a best hit, stopping the run if it is not a
    finite number."""
    try:
        value = float(text)
    except ValueError:
        raise ValueError(f"{path}: coverage '{text}' of the best hit of {query} "
                         "is not a number.") from None
    # float() also accepts 'nan' and 'inf'. Foldseek never prints them, and both
    # would pass the coverage filter of judge_best_hits unnoticed: no comparison
    # with 'nan' is ever true, so it is never "at or below 0.8", and 'inf' is above
    # every cut-off. So they stop the run too, as in the e-value guard.
    if not math.isfinite(value):
        raise ValueError(f"{path}: coverage '{text}' of the best hit of {query} "
                         "is not a finite number.")
    return value


def judge_best_hits(best, rival, catalogue, genes_by_key, hits_path):
    """Apply EcoFoldDB's filter and the borderline flags to every protein's best hit.

    Takes in: the best and rival hits from read_best_hits, the catalogue, and
              {catalogue key: [dictionary genes listing it]}.
    Returns:  {query: structure dict} with
                status    "kept", "failed" or "outside" (best hit outside the
                          catalogue)
                problems  why a catalogue hit failed the filter (empty when kept)
                borderline  the borderline reasons of a kept hit (may be empty)
                genes     the dictionary genes a kept hit gives S evidence to
                target, gene, pathway, key, evalue, qcov, tcov  for the output
              collect_evidence reads it for the S evidence and the structure
              columns of the gene table; main counts the kept hits for the log.
    """
    structure = {}
    for query, hit in best.items():
        entry = catalogue.get(hit.target)
        info = {
            "target": hit.target,
            "gene": entry["gene"] if entry else "",
            "pathway": entry["pathway"] if entry else "",
            "key": entry["key"] if entry else "",
            "evalue": hit.evalue_text,
            "qcov": hit.qcov_text,
            "tcov": hit.tcov_text,
            "problems": [],
            "borderline": [],
            "genes": [],
        }
        structure[query] = info
        # The best hit is a background structure: no S evidence for this protein.
        if entry is None:
            info["status"] = "outside"
            continue

        # EcoFoldDB's filter, with strict inequalities. Every threshold the hit
        # fails is recorded, for the detail of the gene table.
        qcov = _coverage(hit.qcov_text, query, hits_path)
        tcov = _coverage(hit.tcov_text, query, hits_path)
        if hit.evalue >= E_VALUE_MAX:
            info["problems"].append("e-value at or above 1e-10")
        if qcov <= COVERAGE_MIN:
            info["problems"].append("query coverage at or below 0.8")
        if tcov <= COVERAGE_MIN:
            info["problems"].append("target coverage at or below 0.8")
        if info["problems"]:
            info["status"] = "failed"
            continue

        # A kept hit: the genes it gives S evidence to, and its borderline flags.
        info["status"] = "kept"
        info["genes"] = genes_by_key.get(entry["key"], [])
        if hit.evalue >= BORDERLINE_E_VALUE:
            info["borderline"].append("near e-value cut-off")
        if qcov <= BORDERLINE_COVERAGE or tcov <= BORDERLINE_COVERAGE:
            info["borderline"].append("near coverage cut-off")
        # The near-tie test compares the two e-values as Foldseek printed them, in
        # decimal arithmetic. In floating point, 1.000E-23 x 10 comes out just below
        # 1.000E-22, so a rival at exactly ten times the best e-value would miss the
        # flag it should get: the same trap BORDERLINE_E_VALUE avoids for the cut-off.
        other = rival.get(query)
        if other is not None:
            rival_e = decimal.Decimal(other.evalue_text)
            best_e = decimal.Decimal(hit.evalue_text)
            if rival_e <= best_e * BORDERLINE_FACTOR:
                rival_name = other.label if other.label else "a target outside the catalogue"
                info["borderline"].append("near tie with " + rival_name)
    return structure


# ── Reading antiSMASH's regions (supporting evidence only) ───────────────────
# antiSMASH finds biosynthetic gene clusters ("regions") and types them, for
# example NRP-metallophore for a metal-binding molecule, such as a siderophore,
# made by non-ribosomal peptide synthetases. A rule can list a region type as
# supporting evidence; supporting_items reports it when a region of that type
# lies over a protein counted for one of the trait's required genes.

def read_antismash_areas(path, contigs):
    """Read antiSMASH's region list into {contig: [(start, end, products), ...]}.

    Takes in: antiSMASH's {sample}.json (rule secondary_metabolites_analysis)
              and the set of Bakta contig names.
    Does:     converts each area to 1-based inclusive coordinates, as Bakta's
              (antiSMASH writes a 0-based start and an exclusive end). An area
              that runs across the origin of a circular contig keeps a start
              larger than its end: antiSMASH writes the start of the piece at
              the contig end and the end of the piece at the contig start, and
              overlaps() tests the two pieces.
              A record is matched to a Bakta contig by its "name" first and its
              "id" second. antiSMASH reads Bakta's GenBank file, which Bakta
              writes with Biopython. Biopython takes a contig name with one dot
              followed only by digits for an accession with a version number,
              and writes the version as a whole number, dropping leading zeros:
              NODE_1_length_5000_cov_7.05 gets the VERSION line
              NODE_1_length_5000_cov_7.5. antiSMASH reads that line back as the
              record id, while the record name, from the LOCUS line, keeps the
              contig name whole.
    Returns:  (areas by contig, number of areas whose contig Bakta does not have).
              supporting_items reads the areas; main warns about the unmatched.
    """
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    areas = collections.defaultdict(list)
    unmatched = 0
    for record in data.get("records", []):
        # The Bakta contig of this record: by its name first, then by its id.
        contig = None
        for candidate in (record.get("name"), record.get("id")):
            if candidate in contigs:
                contig = candidate
                break
        for area in record.get("areas", []):
            if contig is None:
                unmatched += 1
                continue
            # antiSMASH's 0-based start + 1 is the 1-based first base; its
            # exclusive 0-based end is the 1-based last base, unchanged.
            areas[contig].append((int(area["start"]) + 1, int(area["end"]),
                                  list(area.get("products", []))))
    return areas, unmatched


# ── The two trait tables: reading and checking ───────────────────────────────
# trait_rules_v1.tsv (one row per trait) and trait_genes_v1.tsv (the gene
# dictionary, one row per gene) are read once and checked against each other and
# against the EcoFoldDB catalogue before any genome file is opened, so a slip in
# a table stops the run instead of quietly changing calls. load_trait_tables does
# the checks; the small helpers before it read one cell each.

def read_trait_table(path, columns):
    """Read one trait table into a list of (line number, {column: value}).

    Lines starting with '#' (the provenance line) and blank lines are skipped;
    the first other line is the header. The file has no quoting: a field is
    whatever lies between two tabs. A missing column, or a row with the wrong
    number of fields, stops the run naming the file and the line.
    """
    rows = []
    header = None
    with open(path, "r", encoding="utf-8") as handle:
        for number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                missing = [name for name in columns if name not in header]
                if missing:
                    raise ValueError(f"{path}: line {number}: header is missing "
                                     f"column(s) {', '.join(missing)}.")
                continue
            if len(fields) != len(header):
                raise ValueError(f"{path}: line {number} has {len(fields)} fields, "
                                 f"the header has {len(header)}.")
            rows.append((number, dict(zip(header, fields))))
    if header is None:
        raise ValueError(f"{path}: no header line, so the table is empty.")
    return rows


def split_list(text):
    """'K01428; K01429' -> ['K01428', 'K01429']: the ';'-separated lists of the
    gene dictionary (trait_genes_v1.tsv), with spaces and empty pieces removed."""
    return [item.strip() for item in text.split(";") if item.strip()]


def parse_required(text):
    """'ureA+ureC | ureB+ureC' -> [['ureA', 'ureC'], ['ureB', 'ureC']].

    Alternative sets are separated by '|', genes inside a set by '+'."""
    sets = []
    for alternative in text.split("|"):
        genes = [gene.strip() for gene in alternative.split("+") if gene.strip()]
        if genes:
            sets.append(genes)
    return sets


def parse_colocalised(text):
    """The colocalised column: '' -> None (no condition), 'all' -> 'all' (the whole
    satisfied set), 'a+b; c+d' -> [['a', 'b'], ['c', 'd']] (groups)."""
    text = text.strip()
    if not text:
        return None
    if text.lower() == "all":
        return "all"
    groups = []
    for group in text.split(";"):
        genes = [gene.strip() for gene in group.split("+") if gene.strip()]
        if genes:
            groups.append(genes)
    return groups


def parse_alternatives(text):
    """'speD | nspC' -> ['speD', 'nspC']: the supporting and key_genes columns,
    plain lists of single tokens separated by '|'."""
    return [token.strip() for token in text.split("|") if token.strip()]


def _yes_no(value, where, column, problems):
    """Turn a yes/no cell into True/False, recording a problem for anything else."""
    if value == "yes":
        return True
    if value == "no":
        return False
    problems.append(f"{where}: {column} must be yes or no, not '{value}'")
    return False


def _compile(pattern, where, column, problems):
    """Compile a product pattern case-insensitively; None for an empty cell, and
    None with a problem recorded for a pattern that does not compile."""
    if not pattern:
        return None
    try:
        return re.compile(pattern, re.IGNORECASE)
    except re.error as error:
        problems.append(f"{where}: {column} is not a valid regular expression "
                        f"({error}): {pattern}")
        return None


def load_trait_tables(rules_path, genes_path, catalogue):
    """Read both trait tables and check them against each other and the catalogue.

    Takes in: trait_rules_v1.tsv, trait_genes_v1.tsv and the EcoFoldDB catalogue.
    Checks:   every gene token of required and supporting has exactly one gene
              entry, and every entry is used; no required set names a gene twice
              and no supporting token is listed twice; antiSMASH tokens appear
              only in supporting; colocalised and key_genes genes are in the
              trait's required sets, and each colocalised group lies inside one
              of those sets; every pattern compiles; symbol_needs_product and
              s_alone are yes or no; every ecofolddb key is in the catalogue;
              trait_ids are unique; no required column is empty, and every
              antiSMASH token names a region type.
    Returns:  (rules, genes): rules as a list in table order, genes as
              {gene: entry}. Any failed check stops the run with every problem
              listed, each naming the file, the line and the gene or trait of
              that row (for an entry no trait uses: the file and the entry). A
              broken table would otherwise give a profile of "not called" that
              looks like a result.
    """
    problems = []
    catalogue_keys = {entry["key"] for entry in catalogue.values()}

    # The gene dictionary first, so that every rule can be checked against it.
    genes = {}
    for line, row in read_trait_table(genes_path, GENE_COLUMNS):
        name = row["gene"].strip()
        where = f"{genes_path} line {line} (gene {name})"
        if name in genes:
            problems.append(f"{where}: a second entry for gene {name}")
            continue
        entry = {
            "gene": name,
            "keys": set(split_list(row["ecofolddb"])),
            "kos": set(split_list(row["ko"])),
            "ecs": set(split_list(row["ec"])),
            "symbols": set(split_list(row["bakta_gene"])),
            "symbol_needs_product": _yes_no(row["symbol_needs_product"].strip(),
                                            where, "symbol_needs_product", problems),
            "s_alone": _yes_no(row["s_alone"].strip(), where, "s_alone", problems),
            "include": _compile(row["product_include"], where, "product_include",
                                problems),
            "exclude": _compile(row["product_exclude"], where, "product_exclude",
                                problems),
        }
        for key in sorted(entry["keys"] - catalogue_keys):
            problems.append(f"{where}: ecofolddb key '{key}' is not in the EcoFoldDB "
                            "catalogue (EcoFoldDB_descriptions.txt)")
        genes[name] = entry

    # Then the rules, in table order. used collects every gene a rule names, so
    # that an entry no rule uses can be reported at the end.
    rules = []
    used = set()
    seen_ids = set()
    for line, row in read_trait_table(rules_path, RULE_COLUMNS):
        trait_id = row["trait_id"].strip()
        where = f"{rules_path} line {line} (trait {trait_id})"
        if trait_id in seen_ids:
            problems.append(f"{where}: trait_id used twice")
        seen_ids.add(trait_id)

        sets = parse_required(row["required"])
        if not sets:
            problems.append(f"{where}: required is empty")
        # A gene written twice in one set would be counted twice by the half-of-set
        # test of decide_trait: abcA+abcA+abcB would be "half found" with abcA alone.
        for genes_in_set in sets:
            times_named = collections.Counter(genes_in_set)
            for gene in sorted(times_named):
                if times_named[gene] > 1:
                    problems.append(f"{where}: required set '{'+'.join(genes_in_set)}' "
                                    f"names gene '{gene}' more than once")
        required_genes = {gene for genes_in_set in sets for gene in genes_in_set}
        for gene in sorted(required_genes):
            if gene.startswith(ANTISMASH_PREFIX):
                problems.append(f"{where}: '{gene}' in required; antiSMASH region "
                                "types are supporting evidence only")
            elif gene not in genes:
                problems.append(f"{where}: required gene '{gene}' has no entry in "
                                "the gene dictionary")
            used.add(gene)

        supporting = parse_alternatives(row["supporting"])
        # A token listed twice would be reported twice in the supporting column of
        # the profile; it can only be a slip in the table.
        times_listed = collections.Counter(supporting)
        for token in sorted(times_listed):
            if times_listed[token] > 1:
                problems.append(f"{where}: supporting lists '{token}' more than once")
        for token in supporting:
            if token.startswith(ANTISMASH_PREFIX):
                if not token[len(ANTISMASH_PREFIX):].strip():
                    problems.append(f"{where}: '{token}' names no antiSMASH region type")
                continue
            if token not in genes:
                problems.append(f"{where}: supporting gene '{token}' has no entry in "
                                "the gene dictionary")
            used.add(token)

        colocalised = parse_colocalised(row["colocalised"])
        if isinstance(colocalised, list):
            for group in colocalised:
                outside = [gene for gene in group if gene not in required_genes]
                for gene in outside:
                    problems.append(f"{where}: colocalised gene '{gene}' is not "
                                    "in required")
                # A group applies only to a set that holds all its genes
                # (groups_for_set). A group split over two alternative sets would
                # apply to none of them, and its condition would vanish unseen.
                fits_one_set = any(set(group) <= set(genes_in_set)
                                   for genes_in_set in sets)
                if not outside and not fits_one_set:
                    problems.append(f"{where}: colocalised group '{'+'.join(group)}' "
                                    "is not inside any one required set")

        # key_genes: the genes that set the pathway apart from look-alikes. A
        # partial is shown only if one of them is found somewhere in the genome,
        # counted or as a pseudogene, except a set complete but for pseudogenes,
        # which is shown either way (decide_trait). Each must be one of the
        # trait's required genes.
        key_genes = parse_alternatives(row["key_genes"])
        for gene in key_genes:
            if gene not in required_genes:
                problems.append(f"{where}: key gene '{gene}' (key_genes) is not in "
                                "required")

        rules.append({
            "trait_id": trait_id,
            "trait": row["trait"].strip(),
            "application": row["application"].strip(),
            "caveat": row["caveat"].strip(),
            "sets": sets,
            "supporting": supporting,
            "colocalised": colocalised,
            "key_genes": key_genes,
        })

    for name in genes:
        if name not in used:
            problems.append(f"{genes_path}: gene entry '{name}' is not used by any "
                            "trait")

    if problems:
        raise ValueError("The trait tables did not pass their checks:\n  - "
                         + "\n  - ".join(problems))
    return rules, genes


# ── Evidence: does a protein count for a gene? ───────────────────────────────
# Points 1 to 6 of the module docstring, applied protein by protein in
# collect_evidence with the helpers before it. Out come the gene table and, for
# each dictionary gene, the proteins that count for it and the pseudogenes taken
# as it; the trait decisions below read only those two.

def name_hits(entry, symbol, product, ecs, kos):
    """Say which parts of one annotation name match a dictionary entry.

    Takes in: a gene entry and one protein's symbol, product (without Bakta's
              prefix), EC numbers and KOs (from Bakta, or from Baktfold).
    Does:     the product_include pattern is searched in the product (anywhere,
              ignoring case); a KO or EC number in the entry's lists is a match;
              the exact symbol (case counts) in bakta_gene is a match, except that
              an entry with symbol_needs_product accepts the symbol only together
              with a product or KO match, for symbols Bakta also gives to other
              proteins.
    Returns:  the matching parts, e.g. ["product", "KO K01428"]; empty when the
              name does not point to the gene. The veto is applied by the caller.
    """
    product_hit = bool(entry["include"] and entry["include"].search(product))
    ko_hits = [ko for ko in kos if ko in entry["kos"]]
    ec_hits = [ec for ec in ecs if ec in entry["ecs"]]
    symbol_hit = symbol in entry["symbols"]
    if symbol_hit and entry["symbol_needs_product"] and not (product_hit or ko_hits):
        symbol_hit = False

    parts = []
    if product_hit:
        parts.append("product")
    if symbol_hit:
        parts.append("symbol " + symbol)
    parts += ["KO " + ko for ko in ko_hits]
    parts += ["EC " + ec for ec in ec_hits]
    return parts


def vetoes(entry, product):
    """True when the entry's product_exclude pattern matches this product."""
    return bool(entry["exclude"] and entry["exclude"].search(product))


def structure_note(info, protein_length, max_protein_length, entry_name):
    """Describe, for a name row, what the protein's structure search showed.

    Takes in: the protein's structure dict (None when it had no hit), its length,
              the length limit and the gene of the row.
    Returns:  a short phrase for the detail column of a row where the structure
              gave no S evidence for this gene.
    """
    if info is None:
        if protein_length is not None and protein_length > max_protein_length:
            return (f"not searched for structure: {protein_length} aa, longer than "
                    f"the {max_protein_length} aa limit")
        return "no structure hit"
    if info["status"] == "outside":
        return "best structure hit is outside the catalogue"
    if info["status"] == "failed":
        return "best structure hit failed the filter: " + " and ".join(info["problems"])
    return (f"best structure hit {info['key']} is not in the ecofolddb list of "
            f"{entry_name}")


def protein_row(sample, protein, baktfold, info):
    """Start a gene-table row with the columns that describe the protein itself:
    position, Bakta and Baktfold names, and its best structure hit. The products
    are written as Bakta and Baktfold wrote them, prefix included, so the reader
    sees the mark."""
    row = {
        "sample": sample,
        "locus_tag": protein["locus_tag"],
        "contig": protein["contig"],
        "start": protein["start"],
        "stop": protein["stop"],
        "strand": protein["strand"],
        "cds_index": protein["cds_index"],
        "bakta_gene": protein["gene"],
        "bakta_product": protein["product_written"],
        "ko": ", ".join(protein["kos"]),
        "ec": ", ".join(protein["ecs"]),
        "baktfold_product": baktfold["product_written"] if baktfold else "",
    }
    if info is not None:
        row.update({
            "ecofolddb_target": info["target"],
            "ecofolddb_gene": info["gene"],
            "ecofolddb_pathway": info["pathway"],
            "evalue": info["evalue"],
            "qcov": info["qcov"],
            "tcov": info["tcov"],
        })
    return row


def join_or(items):
    """['rhbC', 'rhbD', 'rhbF'] -> 'rhbC, rhbD or rhbF'; ['speD', 'nspC'] ->
    'speD or nspC'; ['desA'] -> 'desA'. For details that name alternatives."""
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " or " + items[-1]


def evidence_words(item):
    """Say what points to one gene on a protein, naming the gene.

    Takes in: one judged gene row of collect_evidence: its gene, the parts of the
              Bakta and Baktfold names that match it, and whether the kept
              structure hit points to it.
    Returns:  for example 'Bakta name matches pcaH by product, KO K00449', or
              'Bakta name matches abcA by product and the structure hit is to
              abcA'; names first, as in the other details.
    """
    name = item["gene"]
    words = []
    if item["bakta_parts"]:
        words.append(f"Bakta name matches {name} by " + ", ".join(item["bakta_parts"]))
    if item["baktfold_parts"]:
        words.append(f"Baktfold name matches {name} by "
                     + ", ".join(item["baktfold_parts"]))
    if item["has_s"]:
        words.append(f"the structure hit is to {name}")
    return " and ".join(words)


def deciding_source(item):
    """Say which source speaks for a gene when the gene's protein is also named,
    by another source, for another gene of one of its rules (module docstring,
    point 6).

    Takes in: one judged gene row of collect_evidence whose evidence would count:
              its evidence letters (S, Q, F, in that order) and the gene's
              s_alone.
    Returns:  "S" when the structure hit points to the gene and a structure hit
              counts alone for it (s_alone yes); otherwise the gene's first name
              letter: "Q" for Bakta's name, "F" for Baktfold's. A row whose
              evidence would count always has one of the three, because a
              structure hit with no name counts only for a gene with s_alone yes.
    """
    letters = item["letters"]
    if letters[0] == "S" and item["s_alone"]:
        return "S"
    name_letters = [letter for letter in letters if letter != "S"]
    return name_letters[0]


def one_gene_per_rule(judged, traits_of_gene):
    """Find the genes a protein is not given to, because point 6 of the module
    docstring gives it to another gene of one of their rules. A protein is given
    to a gene when it counts for it or, on a protein Bakta marks as a pseudogene,
    when the trait calls take it as a broken copy of that gene.

    Takes in: one protein's judged gene rows from collect_evidence, in the order
              of trait_genes_v1.tsv, each with its gene, whether its evidence
              would count, its evidence letters (S, Q, F, in that order) and the
              gene's s_alone; and {gene: [trait_ids]}, the rules that name each
              gene in their required sets or supporting list.
    Does:     ranks the genes whose evidence would count by the source that
              speaks for each (deciding_source): S, the structure hit, for a gene
              with s_alone yes, then Q, Bakta's name, then F, Baktfold's name;
              genes with the same deciding source keep the order of
              trait_genes_v1.tsv. Going down that list, the protein is given to
              each gene unless the gene shares a rule with a gene it is already
              given to and the two have no letter in common: then no single
              source names both, and the protein stays with the gene that ranks
              first. On a working protein the rows compared are the counted ones;
              on a protein Bakta marks as a pseudogene, the rows with reason
              pseudogene, whose evidence would count but for the mark. A protein
              is either a pseudogene or not, so the two kinds never meet. Rows
              whose evidence would not count anyway take no part:
              structure_only_needs_name rows, and vetoed rows, which
              collect_evidence does not pass in.
    Returns:  {gene the protein is not given to: the judged row of the gene it is
              given to instead}; empty when no two genes clash, as for almost
              every protein. collect_evidence turns the losing row of a working
              protein into counted_for_other_gene, and keeps a pseudogene out of
              the trait calls for the gene it lost.
    """
    source_rank = {"S": 0, "Q": 1, "F": 2}
    in_play = [item for item in judged if item["would_count"]]
    # The ranking: the deciding source of each gene, S before Q before F, so a
    # structure hit to an s_alone gene ranks before Q+F, and Q ranks before the
    # S+F of a gene with s_alone no. Rows with the same deciding source keep the
    # order of trait_genes_v1.tsv, because sorted() leaves tied items in their
    # input order. That order never decides a clash: two rows with the same
    # deciding source share that letter, so neither blocks the other.
    in_play = sorted(in_play, key=lambda item: source_rank[deciding_source(item)])
    still_counting = []        # rows of the genes the protein is given to
    lost_to = {}
    for item in in_play:
        rules_of_item = set(traits_of_gene.get(item["gene"], []))
        winner = None
        for other in still_counting:
            share_a_rule = rules_of_item & set(traits_of_gene.get(other["gene"], []))
            share_a_source = set(item["letters"]) & set(other["letters"])
            if share_a_rule and not share_a_source:
                winner = other
                break
        if winner is None:
            still_counting.append(item)
        else:
            lost_to[item["gene"]] = winner
    return lost_to


def lost_to_detail(loser, winner):
    """The words for a gene row whose protein is given to another gene of one of
    its rules (one_gene_per_rule).

    Takes in: the judged rows of the gene the protein is not given to (the loser)
              and of the gene it is given to (the winner). Their evidence letters
              have none in common, and the winner's deciding source ranks first
              (deciding_source): S, a structure hit to a gene with s_alone yes,
              alone or with a name agreeing, against a name; or Q, a Bakta name,
              against F, a Baktfold name. Either side of Q against F may also
              have a structure hit, to a gene with s_alone no.
    Returns:  what points to the winner and why it ranks first, closed by what that
              means for the protein. On a working protein, whose losing row turns
              counted_for_other_gene, this is the whole detail, opening with what
              points to the loser: "Bakta name matches pcaH by product, KO K00449;
              the structure hit is to pcaG, which counts from structure alone
              (s_alone yes), so the protein counts once, for pcaG". On a
              pseudogene, whose losing row keeps reason pseudogene and a detail
              that already says what points to the loser, it is the note added to
              that detail: "the structure hit is to pcaG, which counts from
              structure alone (s_alone yes), so the broken protein is taken as
              pcaG".
    """
    if winner["letters"] == ["S"]:
        why = ", which counts from structure alone (s_alone yes)"
    elif deciding_source(winner) == "S":
        # The structure hit and a name both point to the winner (S+Q or S+F),
        # and a structure hit counts alone for it.
        why = ", which ranks first as the structure's gene (s_alone yes)"
    elif winner["has_s"]:
        # S+Q for a gene with s_alone no, against F: Bakta's name decides. It is
        # said in a clause of its own, because a "which" right after "the
        # structure hit is to ..." would seem to say that the structure hit ranks
        # first.
        why = "; Bakta's name ranks before a Baktfold name"
    else:
        # Q against F, or against S+F for a gene with s_alone no.
        why = ", which ranks before a Baktfold name"
    if loser["has_s"]:
        # The loser has a structure hit, so its gene has s_alone no (a structure
        # hit to an s_alone gene would have ranked first). Say why it lost anyway.
        why += (f"; the structure hit does not decide, because {loser['gene']} "
                "does not count from structure alone (s_alone no)")
    winner_words = evidence_words(winner) + why
    if loser["reason"] == "pseudogene":
        return f"{winner_words}, so the broken protein is taken as {winner['gene']}"
    return (f"{evidence_words(loser)}; {winner_words}, so the protein counts once, "
            f"for {winner['gene']}")


def collect_evidence(sample, proteins, baktfold_names, structure, lengths, genes,
                     traits_of_gene, max_protein_length):
    """Decide, for every protein and dictionary gene, whether the protein counts.

    Takes in: Bakta's proteins (with positions), the Baktfold names, the judged
              structure hits, the protein lengths, the gene entries, and
              {gene: [trait_ids]}.
    Does:     for each protein, first each gene on its own: the S, Q and F
              evidence; the veto; Bakta's mark; and whether the gene would count
              (see the module docstring). Then the genes whose evidence would
              count are compared (one_gene_per_rule), so that a protein counts for
              two genes of one rule only when one source names both, and a
              pseudogene is taken as two genes of one rule only on the same terms.
              A row is written for every pair with any evidence. A veto row is
              written only when the veto cancelled evidence; an exclude pattern
              that matches an unrelated protein is not news.
    Returns:  (gene-table rows, counted, pseudogenes). counted is
              {gene: [{protein, evidence, borderline}, ...]}, the proteins counted
              for each gene in Bakta's order; borderline holds the reasons when a
              borderline kept hit is part of the evidence. pseudogenes is
              {gene: [{protein, evidence}, ...]}, the proteins taken as broken
              copies of each gene, in Bakta's order: Bakta marks them as
              pseudogenes, their evidence would have counted for the gene (reason
              pseudogene), and one_gene_per_rule gave them to it. decide_trait
              reads it for the pseudogene partial, the key-gene test and the
              details, build_trait_rows for the missing column.
    """
    keyed_rows = []            # (sort key, row), sorted at the end
    counted = collections.defaultdict(list)
    pseudogenes = collections.defaultdict(list)
    for protein in proteins:
        # Rows of this protein sort by contig (in Bakta's order) and CDS index,
        # then by gene, with the rows about the protein itself (gene NA) last.
        position = (protein["contig_order"], protein["cds_index"])
        locus = protein["locus_tag"]
        info = structure.get(locus)
        baktfold = baktfold_names.get(locus)
        length = lengths.get(locus)
        kept = info is not None and info["status"] == "kept"
        structure_genes = set(info["genes"]) if kept else set()
        # The note Bakta's mark adds to the detail of this protein's gene rows;
        # empty for an unmarked protein (BAKTA_MARK_NOTES).
        mark_note = BAKTA_MARK_NOTES.get(protein["mark"], "")

        # The protein's gene rows that are not vetoed. Each is judged on its own in
        # the loop below and finished after it, once one_gene_per_rule has compared
        # the genes whose evidence would count.
        judged = []
        for name, entry in genes.items():
            # Both products are read without Bakta's prefix (read_annotation_table),
            # so a marked protein matches a pattern exactly as an unmarked one does.
            bakta_parts = name_hits(entry, protein["gene"], protein["product"],
                                    protein["ecs"], protein["kos"])
            baktfold_parts = []
            if baktfold is not None:
                baktfold_parts = name_hits(entry, baktfold["gene"], baktfold["product"],
                                           baktfold["ecs"], baktfold["kos"])
            has_s = name in structure_genes
            if not (bakta_parts or baktfold_parts or has_s):
                continue

            # The evidence letters, always in the order S, Q, F.
            letters = []
            if has_s:
                letters.append("S")
            if bakta_parts:
                letters.append("Q")
            if baktfold_parts:
                letters.append("F")
            evidence = "+".join(letters)

            # What each name matched, for the detail column.
            name_notes = []
            if bakta_parts:
                name_notes.append("Bakta name matches by " + ", ".join(bakta_parts))
            if baktfold_parts:
                name_notes.append("Baktfold name matches by " + ", ".join(baktfold_parts))

            row = protein_row(sample, protein, baktfold, info)
            row["gene"] = name
            row["traits"] = ", ".join(traits_of_gene.get(name, []))
            row["evidence"] = evidence
            if has_s:
                row["borderline"] = "; ".join(info["borderline"]) or "no"

            # The veto: the Bakta product first, then the Baktfold product. It cancels
            # a structure hit too, even for an s_alone gene: a fold shared by related
            # enzymes cannot tell them apart, a name can (module docstring, point 3).
            # It comes before Bakta's mark, so a vetoed pseudogene reads
            # vetoed_by_name. The detail quotes the product the pattern was tested
            # on, without Bakta's prefix; the bakta_product column keeps it.
            veto_source = None
            if vetoes(entry, protein["product"]):
                veto_source = ("Bakta", protein["product"])
            elif baktfold is not None and vetoes(entry, baktfold["product"]):
                veto_source = ("Baktfold", baktfold["product"])
            if veto_source is not None:
                detail = (f"product_exclude of {name} matches the {veto_source[0]} "
                          f"product '{veto_source[1]}', so the {evidence} evidence "
                          "does not count")
                if mark_note:
                    detail += "; " + mark_note
                row.update({
                    "counted": "no",
                    "decision": "not counted",
                    "reason": "vetoed_by_name",
                    "detail": detail,
                })
                keyed_rows.append((position + (0, name, row["reason"]), row))
                continue

            # Not vetoed: the evidence gives one of four reasons (module
            # docstring, point 4); only structure_only_needs_name does not count.
            has_name = bool(bakta_parts or baktfold_parts)
            if has_name and has_s:
                reason = "name_and_structure_agree"
                detail = "; ".join(name_notes) + f"; structure hit {info['key']}"
            elif has_name:
                reason = "name_only"
                note = structure_note(info, length, max_protein_length, name)
                detail = "; ".join(name_notes) + "; " + note
            elif entry["s_alone"]:
                reason = "structure_counts_alone"
                detail = (f"structure hit {info['key']}; no name matches, and a "
                          f"structure hit counts alone for {name} (s_alone yes)")
            else:
                reason = "structure_only_needs_name"
                detail = (f"structure hit {info['key']}; no name matches, and {name} "
                          "counts from structure only when a name agrees (s_alone no)")

            would_count = reason != "structure_only_needs_name"
            is_counted = would_count
            # Bakta's mark (module docstring, point 5). A pseudogene whose evidence
            # would count is not counted: its row gets reason pseudogene. Whether
            # the trait calls take it as this gene is settled below, once
            # one_gene_per_rule has compared the protein's genes. A pseudogene that
            # would not count anyway keeps its reason, and every other mark (a gene
            # that runs off an end of its sequence) counts as usual. The mark's
            # note is added to every row when the row is finished below.
            if would_count and protein["mark"] == "pseudogene":
                reason = "pseudogene"
                is_counted = False

            judged.append({
                "gene": name,
                "row": row,
                "reason": reason,
                "detail": detail,
                "would_count": would_count,
                "counted": is_counted,
                "letters": letters,
                "evidence": evidence,
                "bakta_parts": bakta_parts,
                "baktfold_parts": baktfold_parts,
                "has_s": has_s,
                "s_alone": entry["s_alone"],
            })

        # One protein, one gene per rule (module docstring, point 6): which genes
        # this protein is not given to, and the gene it is given to instead. A
        # working protein counts for the genes it is given to; a pseudogene is
        # taken as a broken copy of them.
        lost_to = one_gene_per_rule(judged, traits_of_gene)
        for item in judged:
            name = item["gene"]
            reason = item["reason"]
            detail = item["detail"]
            is_counted = item["counted"]
            if name in lost_to and reason == "pseudogene":
                # A broken protein given to another gene keeps reason pseudogene
                # here (it is broken whichever gene it is); the detail adds the
                # gene it is taken as.
                detail += "; " + lost_to_detail(item, lost_to[name])
            elif name in lost_to:
                reason = "counted_for_other_gene"
                detail = lost_to_detail(item, lost_to[name])
                is_counted = False
            else:
                # A gene the protein is given to names the genes that lost the
                # protein to it, with their evidence, in the order of
                # trait_genes_v1.tsv.
                beaten = []
                for other in judged:
                    other_winner = lost_to.get(other["gene"])
                    if other_winner is not None and other_winner["gene"] == name:
                        beaten.append(f"{other['gene']} ({other['evidence']})")
                if beaten and reason == "pseudogene":
                    detail += (f"; the broken protein is taken as {name}, not as "
                               + join_or(beaten))
                elif beaten:
                    detail += (f"; the protein counts once, for {name}, not for "
                               + join_or(beaten))
                # A pseudogene given to this gene is the broken copy the trait
                # calls see: decide_trait, broken_outside and the profile's missing
                # column read pseudogenes. A gene that lost the broken protein to
                # another gene (the first branch above) gets no entry for it in
                # pseudogenes.
                if reason == "pseudogene":
                    pseudogenes[name].append({"protein": protein,
                                              "evidence": item["evidence"]})
            # The note of Bakta's mark ends the detail; on a pseudogene row it also
            # says why the row does not count.
            if reason == "pseudogene":
                detail += "; " + mark_note + ", so it does not count"
            elif mark_note:
                detail += "; " + mark_note

            row = item["row"]
            row.update({
                "counted": "yes" if is_counted else "no",
                "decision": "counted" if is_counted else "not counted",
                "reason": reason,
                "detail": detail,
            })
            keyed_rows.append((position + (0, name, reason), row))
            # The proteins counted for each gene, in Bakta's order, with the
            # borderline reasons of a structure hit behind the count.
            if is_counted:
                counted[name].append({
                    "protein": protein,
                    "evidence": item["evidence"],
                    "borderline": info["borderline"] if item["has_s"] else [],
                })

        # Rows about the protein rather than one gene (gene NA).
        if kept and not info["genes"]:
            row = protein_row(sample, protein, baktfold, info)
            row.update({
                "counted": "no",
                "evidence": "S",
                "borderline": "; ".join(info["borderline"]) or "no",
                "decision": "not counted",
                "reason": "catalogue_gene_unused",
                "detail": (f"kept structure hit {info['key']}, which no gene of "
                           "the gene dictionary uses"),
            })
            keyed_rows.append((position + (1, "", row["reason"]), row))
        # A protein over the length limit was left out of the structure search
        # (rule trait_structure_query), so it cannot have hits; the row says so.
        # Should a hit table from another set-up hold hits for it anyway, it was
        # searched, and the row would be wrong, hence "info is None".
        if info is None and length is not None and length > max_protein_length:
            row = protein_row(sample, protein, baktfold, info)
            row.update({
                "counted": "no",
                "decision": "not searched",
                "reason": "not_searched_too_long",
                "detail": (f"{length} aa, longer than the {max_protein_length} aa "
                           "limit of the structure search; its names were still "
                           "checked"),
            })
            keyed_rows.append((position + (1, "", row["reason"]), row))

    keyed_rows.sort(key=lambda pair: pair[0])
    rows = [row for _key, row in keyed_rows]
    return rows, counted, pseudogenes


# ── Deciding each trait ──────────────────────────────────────────────────────
# "How a trait is decided" in the module docstring. decide_trait applies one rule
# to the counted proteins (chained and groups_for_set test co-localisation);
# describe_gene, describe_missing and broken_outside write genes the way the
# profile and the audit show them; supporting_items gathers the supporting
# evidence (overlaps places antiSMASH regions); build_trait_rows makes one profile
# row and one audit row per trait.

def chained(counted, genes, max_gap):
    """True when one chain on one contig holds a counted protein of every gene.

    Takes in: the counted proteins per gene, the genes of one co-localised group
              and the largest allowed gap in CDS steps.
    Does:     per contig, walks the counted proteins of the group's genes in CDS
              order and starts a new chain whenever the step from the previous
              protein is larger than max_gap. Each gene must therefore lie at most
              max_gap genes from the one before it, not from the first.
    Returns:  True as soon as a chain holds every gene of the group.
    """
    wanted = set(genes)
    # (CDS index, gene) of every counted protein of the group's genes, by contig.
    items_by_contig = collections.defaultdict(list)
    for gene in genes:
        for found in counted.get(gene, []):
            protein = found["protein"]
            items_by_contig[protein["contig"]].append((protein["cds_index"], gene))
    for items in items_by_contig.values():
        items.sort()
        chain = set()
        previous_index = None
        for index, gene in items:
            if previous_index is not None and index - previous_index > max_gap:
                chain = set()
            chain.add(gene)
            previous_index = index
            if wanted <= chain:
                return True
    return False


def groups_for_set(colocalised, genes_in_set):
    """The co-localisation groups that apply to one required set.

    "all" means the whole set; a listed group applies only when all its genes are
    in this set, so a group written for another alternative set is ignored here.
    """
    if colocalised is None:
        return []
    if colocalised == "all":
        return [list(genes_in_set)]
    return [group for group in colocalised if set(group) <= set(genes_in_set)]


def decide_trait(rule, counted, pseudogenes, max_gap):
    """Apply one trait's rule to the genes counted in this genome.

    Takes in: one rule from load_trait_tables; the counted proteins per gene and
              the pseudogenes taken as each gene (both from collect_evidence); the
              largest co-localisation gap.
    Does:     the decision of the module docstring, set by set in table order:
              the first complete, co-localised set is the call; failing that, the
              first set that qualifies as a partial and is not withheld is the
              partial. A set qualifies when it has two or more genes, all counted
              but not co-localised; or when each of its genes is counted or a
              pseudogene, at least one a pseudogene (a set of one gene too); or
              when it has two or more genes, at least half of them counted.
              key_genes can withhold the first and the last kind, never the
              middle one. A called set wins over any partial, wherever it stands
              in the table. When every partial is withheld, the trait is not
              called and the detail shows the first one. The detail of a partial
              or of a trait not called also names the trait's genes outside the
              set it reports that are found only as a pseudogene (broken_outside).
    Returns:  {decision, reason, detail, call_set}; call_set is the satisfied or
              partial set, None when not called.
    """
    # present: the genes with at least one counted protein. broken: the genes with
    # at least one pseudogene taken as them (collect_evidence).
    present = {gene for gene, found in counted.items() if found}
    broken = {gene for gene, found in pseudogenes.items() if found}
    # The rule's key genes are found when one of them is counted anywhere in the
    # genome, or found there as a pseudogene (a pseudogene is taken as it): the
    # gene is then in the genome, even if broken. A rule without key_genes has
    # nothing to find.
    key_gene_found = (not rule["key_genes"]
                      or any(gene in present or gene in broken
                             for gene in rule["key_genes"]))
    # The words that open the note naming broken genes of the other sets.
    outside_note = "; outside this set, found only as a pseudogene: "
    partial = None      # the first set that qualifies as a partial and stands
    withheld = None     # the first partial held back for want of a key gene
    for genes_in_set in rule["sets"]:
        set_text = "+".join(genes_in_set)
        found = [gene for gene in genes_in_set if gene in present]
        missing = [gene for gene in genes_in_set if gene not in present]
        groups = groups_for_set(rule["colocalised"], genes_in_set)

        failed_groups = []
        if not missing:
            failed_groups = [group for group in groups
                             if not chained(counted, group, max_gap)]
            if not failed_groups:
                detail = set_text + " complete"
                if groups:
                    detail += "; co-localised: " + "; ".join("+".join(group)
                                                            for group in groups)
                return {"decision": "called", "reason": "set_complete",
                        "detail": detail, "call_set": genes_in_set}

        # Only the first partial that stands is kept. The loop still goes on,
        # because a later set may be called.
        if partial is not None:
            continue
        # The genes of this set that are not counted, as the detail writes them:
        # gene:locus_tag(pseudogene) for a gene found only as a pseudogene, the
        # bare name for any other (describe_missing).
        missing_text = ", ".join(describe_missing(gene, pseudogenes)
                                 for gene in missing)
        missing_broken = [gene for gene in missing if gene in broken]
        # The three partials, tested in this order. A set of one gene can be
        # neither half found nor split over the genome, so the first and the last
        # test need two or more genes.
        candidate = None
        if not missing and len(genes_in_set) >= 2:
            failed_text = "; ".join("+".join(group) for group in failed_groups)
            candidate = {"decision": "partial",
                         "reason": "set_complete_not_colocalised",
                         "detail": (f"{set_text} complete but not co-localised "
                                    f"({failed_text}: not within {max_gap} genes of "
                                    "each other on one contig)"),
                         "call_set": genes_in_set}
        elif missing and len(missing_broken) == len(missing):
            # Every gene not counted is a pseudogene: the whole set is in the
            # genome, but part of it is broken. This holds for a set of one gene
            # too, whose gene is there but broken. Co-localisation is not tested
            # here: the set is partial wherever its genes lie.
            candidate = {"decision": "partial",
                         "reason": "set_complete_but_pseudogene",
                         "detail": f"{set_text} complete except {missing_text}",
                         "call_set": genes_in_set}
        elif len(genes_in_set) >= 2 and len(found) >= math.ceil(len(genes_in_set) / 2):
            candidate = {"decision": "partial",
                         "reason": "half_of_set_found",
                         "detail": (f"missing {missing_text} ({len(found)} of "
                                    f"{len(genes_in_set)} genes of {set_text} "
                                    "found)"),
                         "call_set": genes_in_set}
        if candidate is None:
            continue

        # key_genes is there to hold back partials that a related pathway,
        # sharing some of the genes, could produce: an aerobactin cluster also
        # matches desB, desC and desD of desferrioxamine, but holds no desA, the
        # key gene. It never holds back a set that is complete except for
        # pseudogenes: every gene of that set is in the genome, and the partial is
        # there to show the broken one. A set of one gene makes the point: when
        # the rule's key genes all belong to its other routes, a broken gene of
        # that one-gene route would otherwise read "not called". A withheld
        # partial is kept aside, and the loop goes on, since a later set may
        # still stand.
        if key_gene_found or candidate["reason"] == "set_complete_but_pseudogene":
            partial = candidate
        elif withheld is None:
            withheld = candidate

    # A partial stands: its detail also names the broken genes outside its set.
    if partial is not None:
        others = broken_outside(rule, present, pseudogenes, partial["call_set"])
        if others:
            partial["detail"] += outside_note + ", ".join(others)
        return partial
    # Every partial was withheld: not called, and the detail shows the first one.
    if withheld is not None:
        withheld_detail = withheld["detail"]
        others = broken_outside(rule, present, pseudogenes, withheld["call_set"])
        if others:
            withheld_detail += outside_note + ", ".join(others)
        # "no key gene (desA) found (...)", or with several key genes
        # "no key gene (rhbC, rhbD or rhbF) found (...)".
        return {"decision": "not called", "reason": "partial_withheld",
                "detail": (f"no key gene ({join_or(rule['key_genes'])}) found "
                           f"({withheld_detail})"),
                "call_set": None}

    # Not called: report the closest set, the one with most genes found (the
    # first one on a tie), so the reader sees how far the genome is from a call.
    closest = None
    closest_found = []
    for genes_in_set in rule["sets"]:
        found = [gene for gene in genes_in_set if gene in present]
        if closest is None or len(found) > len(closest_found):
            closest = genes_in_set
            closest_found = found
    if not closest_found:
        # No gene of any set counts. A required gene found only as a pseudogene
        # is then the one trace of the trait in this genome, so the detail names
        # it. That detail says "counted": "no gene found" beside "found only as a
        # pseudogene" would read as a contradiction.
        others = broken_outside(rule, present, pseudogenes, [])
        if others:
            detail = ("no gene of any required set counted; found only as a "
                      "pseudogene: " + ", ".join(others))
        else:
            detail = "no gene of any required set found"
        return {"decision": "not called", "reason": "no_gene_found",
                "detail": detail, "call_set": None}
    missing = [gene for gene in closest if gene not in present]
    missing_text = ", ".join(describe_missing(gene, pseudogenes) for gene in missing)
    detail = (f"{len(closest_found)} of {len(closest)} genes found in the closest "
              f"set {'+'.join(closest)} (missing {missing_text})")
    others = broken_outside(rule, present, pseudogenes, closest)
    if others:
        detail += outside_note + ", ".join(others)
    return {"decision": "not called", "reason": "too_few_genes",
            "detail": detail, "call_set": None}


def describe_gene(gene, found):
    """'ureC:SAMPLE1_00012(S+Q),SAMPLE1_00340(Q)': a gene with its counted proteins."""
    proteins = ",".join(f"{item['protein']['locus_tag']}({item['evidence']})"
                        for item in found)
    return f"{gene}:{proteins}"


def describe_missing(gene, pseudogenes):
    """A gene that is not counted, as the audit detail and the profile's missing
    column write it.

    Takes in: the gene and the pseudogenes taken as each gene, from
              collect_evidence.
    Returns:  'pvdA:SAMPLE1_01234(pseudogene)' when a pseudogene is taken as the
              gene (several proteins joined by ',', as describe_gene joins them),
              or the bare name 'pvdA' otherwise: when no protein matched the
              gene, or none that counts for it or is taken as it (a vetoed one,
              say, or a pseudogene given to another gene).
    """
    broken_copies = pseudogenes.get(gene)
    if not broken_copies:
        return gene
    loci = ",".join(f"{item['protein']['locus_tag']}(pseudogene)"
                    for item in broken_copies)
    return f"{gene}:{loci}"


def broken_outside(rule, present, pseudogenes, genes_named):
    """List the trait's required genes found only as a pseudogene outside one set.

    Takes in: the rule; the genes counted in this genome (present); the
              pseudogenes taken as each gene, from collect_evidence; the genes the
              detail already speaks for: the set it reports, whose missing genes
              it names itself, or an empty list when it reports no set.
    Does:     walks the rule's required genes in table order and keeps each one
              once when it is not in genes_named, no protein counts for it, and a
              pseudogene is taken as it. A gene with a working copy is left out:
              it counts, so nothing about it is hidden.
    Returns:  those genes as describe_missing writes them, for example
              ['abcD:SAMPLE1_00004(pseudogene)']; an empty list when there are
              none.
    Why:      the detail of a decision names only the genes of the set it reports.
              Without this list, a gene broken in another set, or the one broken
              gene of a trait whose other genes are all absent, would show in the
              gene table and nowhere in the audit.
    """
    texts = []
    seen = set(genes_named)
    for genes_in_set in rule["sets"]:
        for gene in genes_in_set:
            if gene in seen:
                continue
            seen.add(gene)
            if gene not in present and pseudogenes.get(gene):
                texts.append(describe_missing(gene, pseudogenes))
    return texts


def overlaps(protein, area_start, area_end):
    """True when a protein and an antiSMASH area share at least one base.

    A protein whose start is larger than its stop runs across the end of a
    circular contig and covers start..contig end plus 1..stop; it overlaps an area
    that reaches either piece. An area can run across the end in the same way
    (start larger than end); it is tested as its two pieces. The contig length is
    not needed for that: no protein lies past the contig end, so the piece at the
    contig end can simply run to infinity.
    """
    if area_start > area_end:
        at_contig_end = overlaps(protein, area_start, math.inf)
        at_contig_start = overlaps(protein, 1, area_end)
        return at_contig_end or at_contig_start
    start = protein["start"]
    stop = protein["stop"]
    if start <= stop:
        return start <= area_end and stop >= area_start
    return area_end >= start or area_start <= stop


def supporting_items(rule, counted, areas):
    """List the supporting evidence found for one trait.

    Takes in: the rule, the counted proteins per gene, and the antiSMASH areas by
              contig (None when no antiSMASH output was given).
    Does:     reports each counted supporting gene with its proteins, and each
              antiSMASH region type of the supporting list that antiSMASH called
              over a protein counted for one of the trait's required genes, on the
              same contig. A region type must match exactly ("NRP-metallophore" is
              not "NRPS"). Without antiSMASH output the type is "not checked".
    Returns:  a list of strings, in the order of the supporting column.
    """
    # The proteins counted for the trait's required genes, by contig: an antiSMASH
    # region is reported only when it lies over one of them.
    required_genes = {gene for genes_in_set in rule["sets"] for gene in genes_in_set}
    required_proteins = collections.defaultdict(list)
    for gene in sorted(required_genes):
        for found in counted.get(gene, []):
            protein = found["protein"]
            required_proteins[protein["contig"]].append(protein)

    items = []
    for token in rule["supporting"]:
        # A supporting gene: listed with its counted proteins, when it has any.
        if not token.startswith(ANTISMASH_PREFIX):
            if counted.get(token):
                items.append(describe_gene(token, counted[token]))
            continue
        if areas is None:
            items.append(f"{token}(not checked)")
            continue
        # An antiSMASH region type: every region of exactly that type that lies
        # over a protein counted for a required gene, on the same contig.
        region_type = token[len(ANTISMASH_PREFIX):]
        places = []
        for contig, contig_areas in areas.items():
            for area_start, area_end, products in contig_areas:
                if region_type not in products:
                    continue
                if any(overlaps(protein, area_start, area_end)
                       for protein in required_proteins.get(contig, [])):
                    places.append(f"{contig}:{area_start}-{area_end}")
        if places:
            items.append(f"{token}({','.join(places)})")
    return items


def build_trait_rows(sample, rules, counted, pseudogenes, areas, max_gap):
    """Decide every trait and build its profile row and its audit row.

    Takes in: the rules in table order, the counted proteins per gene and the
              pseudogenes per gene (both from collect_evidence), the antiSMASH
              areas (or None) and the co-localisation gap.
    Does:     the genes, borderline, structure_only and supporting columns are
              built from counted proteins only; a pseudogene shows only in the
              missing column (and in the audit detail).
    Returns:  (profile rows, audit rows), one of each per trait, in table order.
    """
    profile_rows = []
    audit_rows = []
    for rule in rules:
        verdict = decide_trait(rule, counted, pseudogenes, max_gap)
        row = {
            "sample": sample,
            "trait_id": rule["trait_id"],
            "trait": rule["trait"],
            "application": rule["application"],
            "decision": verdict["decision"],
            "supporting": "; ".join(supporting_items(rule, counted, areas)),
            "caveat": rule["caveat"],
        }
        # The genes, missing, borderline and structure_only columns describe the
        # called or partial set only; a trait not called leaves them NA.
        call_set = verdict["call_set"]
        if call_set is not None:
            call_genes = [gene for gene in call_set if counted.get(gene)]
            row["genes"] = "; ".join(describe_gene(gene, counted[gene])
                                     for gene in call_genes)
            row["missing"] = ", ".join(describe_missing(gene, pseudogenes)
                                       for gene in call_set if not counted.get(gene))
            any_borderline = any(found["borderline"]
                                 for gene in call_genes for found in counted[gene])
            row["borderline"] = "yes" if any_borderline else "no"
            row["structure_only"] = ", ".join(
                gene for gene in call_genes
                if all(found["evidence"] == "S" for found in counted[gene]))
        profile_rows.append(row)
        audit_rows.append({
            "sample": sample,
            "trait_id": rule["trait_id"],
            "decision": verdict["decision"],
            "reason": verdict["reason"],
            "detail": verdict["detail"],
        })
    return profile_rows, audit_rows


def traits_using_each_gene(rules):
    """{gene: [trait_ids]}: the traits whose required or supporting genes name it,
    in table order. collect_evidence writes it in the traits column of the gene
    table, and one_gene_per_rule uses it to compare only genes that share a
    rule."""
    traits = collections.defaultdict(list)
    for rule in rules:
        genes = [gene for genes_in_set in rule["sets"] for gene in genes_in_set]
        genes += [token for token in rule["supporting"]
                  if not token.startswith(ANTISMASH_PREFIX)]
        for gene in genes:
            if rule["trait_id"] not in traits[gene]:
                traits[gene].append(rule["trait_id"])
    return traits


# ── Writing the three tables and the command line ────────────────────────────
# main() reads the catalogue and checks the two tables first, then reads the
# genome, the structure hits and antiSMASH's regions, runs the steps above and
# writes the three tables.

def _cell(value):
    """One output cell: NA for nothing, and no tab or line break inside a cell, so
    every row stays one line whichever way the table is read."""
    if value is None:
        return "NA"
    text = str(value)
    if text == "":
        return "NA"
    return text.replace("\t", " ").replace("\r", " ").replace("\n", " ")


def write_tsv(path, columns, rows):
    """Write rows (dicts) as a tab-separated table with a header.

    The header is written even when there are no rows: an empty table with a
    header says "we looked and found nothing"; a missing file says nothing.
    No cell is quoted: _cell keeps tabs and line breaks out, so a plain split on
    tabs reads every row, in a spreadsheet, awk or Python alike.
    """
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(_cell(row.get(column)) for column in columns) + "\n")


def file_fingerprint(path):
    """Fingerprint one file: {"file": its name, "sha256": its SHA-256 digest}.

    The SHA-256 digest is a 64-character code computed from every byte of the
    file; any change to the file, even one character, gives a different code.
    Only the file name is kept beside it, as a label for a person reading the
    JSON: the folder would reveal where the analysis ran, and comparing
    contents rather than paths keeps the record valid when the output folder
    is moved.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"file": os.path.basename(path), "sha256": digest.hexdigest()}


def write_provenance(path, sample, inputs, parameters, outputs):
    """Write {sample}_trait_provenance.json: what these calls were made with.

    Takes in: the output path; the sample name; the fingerprints of the three
              catalogue files taken when main() started (inputs); the calling
              parameters; and the paths of the three tables just written
              (outputs), which are fingerprinted here.
    Why: the catalogue's caveats, assays and references are quoted beside the
    calls by rule traits_summary (in the MultiQC report) and rule
    traits_workbook (in the traits workbook). Both compare these fingerprints
    with the files they read, so a catalogue edited after the calls were made
    is caught instead of being quoted beside calls it did not make.
    Written last, after the three tables, so a run that stops half way leaves
    no record vouching for incomplete tables.
    """
    record = {
        "schema_version": "1.0",
        "sample": sample,
        "producer": file_fingerprint(__file__),
        "inputs": inputs,
        "parameters": parameters,
        "outputs": {name: file_fingerprint(table) for name, table in outputs.items()},
    }
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(record, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")


def main(argv=None):
    """Command-line entry point used by rule trait_profile.

    All three tables are always written, so the rule's outputs exist even for a
    genome where nothing is found.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Call a strain's trait profile from sets of genes, using EcoFoldDB "
            "structure hits, Bakta names and, when given, Baktfold names, with a "
            "reason for every decision."
        )
    )
    parser.add_argument("--sample", required=True,
                        help="Sample name; the first column of every output.")
    parser.add_argument("--bakta-tsv", required=True,
                        help="Bakta's {sample}.tsv annotation table.")
    parser.add_argument("--bakta-faa", required=True,
                        help="Bakta's {sample}.faa, read for protein lengths.")
    parser.add_argument("--hits", required=True,
                        help="Foldseek hits against EcoFoldDB (query, target, "
                             "evalue, qcov, tcov; no header); plain or .gz.")
    parser.add_argument("--catalogue", required=True,
                        help="EcoFoldDB_descriptions.txt of the EcoFoldDB database.")
    parser.add_argument("--rules", default=DEFAULT_RULES_TSV,
                        help="Trait rules table. Default: trait_rules_v1.tsv beside "
                             "this script.")
    parser.add_argument("--genes", default=DEFAULT_GENES_TSV,
                        help="Trait gene dictionary. Default: trait_genes_v1.tsv "
                             "beside this script.")
    parser.add_argument("--antismash-json", default=None,
                        help="antiSMASH {sample}.json (optional; supporting "
                             "evidence only).")
    parser.add_argument("--baktfold-tsv", default=None,
                        help="Baktfold {sample}.tsv (optional; only when Baktfold "
                             "ran).")
    parser.add_argument("--max-gene-gap", type=int, default=DEFAULT_MAX_GENE_GAP,
                        help="Co-localised genes: each at most this many genes "
                             "from the previous one; 1 or more. Default "
                             f"{DEFAULT_MAX_GENE_GAP}.")
    parser.add_argument("--max-protein-length", type=int,
                        default=DEFAULT_MAX_PROTEIN_LENGTH,
                        help="Proteins longer than this were not searched for "
                             f"structure. Default {DEFAULT_MAX_PROTEIN_LENGTH}.")
    parser.add_argument("--out-genes", required=True,
                        help="Output: one row per protein and gene with evidence.")
    parser.add_argument("--out-profile", required=True,
                        help="Output: one row per trait.")
    parser.add_argument("--out-audit", required=True,
                        help="Output: the reason behind each trait decision.")
    parser.add_argument("--out-provenance",
                        help="Optional: write the fingerprints of the catalogue "
                             "files and output tables, with the calling "
                             "parameters, to this JSON file.")
    args = parser.parse_args(argv)
    # Next-door genes are 1 apart, so with a gap of 0 no two different proteins
    # could ever be co-localised, and every trait that asks for it would be
    # partial at best. traits.max_gene_gap reaches the script only through this
    # option.
    if args.max_gene_gap < 1:
        parser.error("--max-gene-gap must be 1 or more")
    if args.max_protein_length < 1:
        parser.error("--max-protein-length must be 1 or more")

    # Fingerprint the three catalogue files before reading them, so the record
    # describes the files these calls were made from. A run by hand may leave
    # out --out-provenance; the workflow always asks for it.
    provenance_inputs = None
    if args.out_provenance:
        provenance_inputs = {
            "rules": file_fingerprint(args.rules),
            "gene_dictionary": file_fingerprint(args.genes),
            "catalogue": file_fingerprint(args.catalogue),
        }

    # 1. The two tables, checked against each other and the catalogue before
    #    anything else is read.
    catalogue = read_catalogue(args.catalogue)
    rules, genes = load_trait_tables(args.rules, args.genes, catalogue)
    # Which dictionary genes list each catalogue key: a kept hit to that entry
    # gives S evidence to all of them (judge_best_hits).
    genes_by_key = collections.defaultdict(list)
    for name, entry in genes.items():
        for key in entry["keys"]:
            genes_by_key[key].append(name)

    # 2. The genome: Bakta's proteins with positions and lengths, Baktfold names.
    proteins = read_annotation_table(args.bakta_tsv)
    add_positions(proteins)
    lengths = read_protein_lengths(args.bakta_faa)
    baktfold_names = {}
    if args.baktfold_tsv:
        baktfold_names = read_baktfold_names(args.baktfold_tsv, proteins)

    # 3. The structure evidence. A hit table whose proteins are not this genome's
    #    means the wrong file was passed: without the check, every trait would
    #    quietly be called from names alone.
    best, rival = read_best_hits(args.hits, catalogue)
    known = {protein["locus_tag"] for protein in proteins}
    unknown = [query for query in best if query not in known]
    if best and len(unknown) == len(best):
        raise ValueError(f"{args.hits}: none of its {len(best)} query proteins is in "
                         f"{args.bakta_tsv}, so the hits belong to another genome.")
    if unknown:
        print(f"WARNING: {len(unknown)} protein(s) in the hit table are not in the "
              f"Bakta table and were ignored, e.g. {unknown[0]}.", file=sys.stderr)
    # An empty hit table is accepted, so that a hand-run without the structure
    # search still works, but then every trait is decided from names alone. The
    # search of a whole genome against EcoFoldDB gives about 3 million hits, so an
    # empty table for a genome with proteins short enough to be searched most
    # likely means the search failed or the wrong file was given. A warning, not a
    # stop.
    if not best:
        searchable = 0
        for protein in proteins:
            length = lengths.get(protein["locus_tag"])
            if length is not None and length <= args.max_protein_length:
                searchable += 1
        if searchable:
            print(f"WARNING: the hit table {args.hits} holds no hits, although "
                  f"{searchable} protein(s) are short enough for the structure "
                  f"search (at most {args.max_protein_length} aa). Every trait is "
                  "decided from names alone.", file=sys.stderr)
    structure = judge_best_hits(best, rival, catalogue, genes_by_key, args.hits)

    # antiSMASH's regions, when given. None makes every antiSMASH supporting type
    # read "(not checked)" (supporting_items).
    areas = None
    if args.antismash_json:
        contigs = {protein["contig"] for protein in proteins}
        areas, unmatched = read_antismash_areas(args.antismash_json, contigs)
        if unmatched:
            print(f"WARNING: {unmatched} antiSMASH region(s) lie on contigs the "
                  "Bakta table does not name and were not compared.",
                  file=sys.stderr)

    # 4. Evidence per protein and gene, then one decision per trait.
    traits_of_gene = traits_using_each_gene(rules)
    gene_rows, counted, pseudogenes = collect_evidence(
        args.sample, proteins, baktfold_names, structure, lengths, genes,
        traits_of_gene, args.max_protein_length)
    profile_rows, audit_rows = build_trait_rows(args.sample, rules, counted,
                                                pseudogenes, areas, args.max_gene_gap)

    # 5. The three tables, written even when they hold no rows.
    write_tsv(args.out_genes, GENE_TABLE_COLUMNS, gene_rows)
    write_tsv(args.out_profile, PROFILE_COLUMNS, profile_rows)
    write_tsv(args.out_audit, AUDIT_COLUMNS, audit_rows)

    # The provenance record last, once the three tables exist (see
    # write_provenance for why).
    if args.out_provenance:
        write_provenance(args.out_provenance, args.sample, provenance_inputs, {
            "max_gene_gap": args.max_gene_gap,
            "max_protein_length": args.max_protein_length,
            "structure_evalue_max": E_VALUE_MAX,
            "structure_coverage_min": COVERAGE_MIN,
            "borderline_evalue": BORDERLINE_E_VALUE,
            "borderline_coverage": BORDERLINE_COVERAGE,
            "borderline_factor": BORDERLINE_FACTOR,
            "antismash_regions_used": args.antismash_json is not None,
            "baktfold_names_used": args.baktfold_tsv is not None,
        }, {
            "genes": args.out_genes,
            "profile": args.out_profile,
            "audit": args.out_audit,
        })

    # Closing summary. The rule sends stdout to its log, so the counts and the
    # path to the reasons are there even when nothing is called.
    decisions = collections.Counter(row["decision"] for row in profile_rows)
    n_kept = sum(1 for info in structure.values() if info["status"] == "kept")
    print(f"Sample {args.sample}: {decisions['called']} trait(s) called, "
          f"{decisions['partial']} partial, {decisions['not called']} not called, "
          f"of {len(rules)}. Structure: {len(structure)} protein(s) with hits, "
          f"{n_kept} kept by the EcoFoldDB filter. Baktfold names: "
          f"{len(baktfold_names)}. Reasons in {args.out_audit} and {args.out_genes}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
