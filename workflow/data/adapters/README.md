# Adapter sequences for ONT trimming and the adapter screen

Two steps read the files in this folder:

- rule `find_ont_adapter` (`workflow/rules/shared/05_ont_adapters.smk`) names the ONT
  barcodes it finds in the reads from `ont_chtd.fasta`, and trims every barcode of a file
  that holds several;
- rule `adapter_search_input` (`workflow/rules/shared/20_qc.smk`) joins both files into the
  list of adapters the finished genome is searched for. Two more sets are added there that
  are not shipped: the Illumina records of NCBI UniVec, which rule `univec_illumina`
  downloads into the run (UniVec notes on every Illumina record that the oligonucleotide
  sequence is copyright Illumina, Inc.), and, for an ONT sample whose reads were trimmed,
  the sequences that were trimmed.

| File | Records | Source |
|---|---|---|
| `ont_chtd.fasta` | 196 | Oxford Nanopore Chemistry Technical Document, version CHTD_500_v1_revAV_24July2026 (<https://nanoporetech.com/document/chemistry-technical-document>), sections 14 and 15: the adapter top strand (shared by the kit 14 ligation, native and rapid adapters), the ligation and native adapter bottom strands, the 3' flanking sequence of the rapid barcodes, the 96 native barcodes NB01-96 between their flanks (AAGGTTAA-barcode-CAGCACCT) and the 96 rapid barcodes RB01-96 between theirs (GCTTGGGTGTTTAACC-barcode-3' flank). ONT publishes the rapid barcodes only in part. |
| `univec_pacbio_mgi.fasta` | 5 | NCBI UniVec, the file of 2024-08-21 at <https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec>: the PacBio blunt (SMRTbell) adapter NGB00972.1, the PacBio C2 primer NGB00973.1, the PacBio ULI amplification adapter NGB03000.1 and the two MGI DNBSEQ adaptors NGB05000.1 and NGB05001.1. UniVec is a work of the US Government. |

Every FASTA header names its source. The native barcodes were read from the document's
table, and its forward and reverse columns were checked to be reverse complements of each
other. All 192 barcode sequences also appear, on one strand or the other, in the adapter list
of Porechop (<https://github.com/rrwick/Porechop>, `porechop/adapters.py`), an independent
copy of ONT's published barcodes. Three of the five UniVec records (NGB03000.1, NGB05000.1,
NGB05001.1) are newer than UniVec build 10.1, the build rule `univec_illumina` downloads, which
is one reason all five are shipped rather than taken from that download.

One sequence is not synthetic. The 3' flank of the rapid barcodes, the 50 bp after every
rapid barcode (`ONT_rapid_barcode_3prime_flank`, and the last 50 bp of every RB record), is
identical to the last 50 bp of the phage Mu genome (NC_000929.1), Mu's right end, which a Mu
transposase joins to DNA. A bacterial genome that carries a Mu prophage holds the same 50 bp,
so the adapter screen lists a match to it alone inside a contig but does not count it
(`workflow/scripts/20_qc/adapter_screen.py`), and the read scan finds a rapid barcode through
it only when a published barcode sits in front of it
(`workflow/scripts/05_ont_adapters/find_ont_adapter.py`).

Both files are read through rule `params`, not as inputs: a file in the repository gets a new
modification time at every checkout, and as an input that alone would re-run the trimming,
and with it the assembly, of every long-read sample in an existing output folder. A change to
these files therefore does not re-run anything by itself.
