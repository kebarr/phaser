# Phaser sample test dataset

A small synthetic GFA + 10x-style mapping file, built to exercise the
whole phasing pipeline (bubble finding, haplotype enumeration, barcode
scoring, and the noise-filtering heuristics) end to end.

## Files

- `sample.gfa` — a linear assembly graph with two diploid bubbles:
  `homA -> {het1a, het1b} -> homB -> {het2a, het2b} -> homC`.
  `homA`/`homB`/`homC` are the homozygous backbone; each bubble
  represents one heterozygous site with two allelic contigs.
- `graphs.txt` — the graph-list input: one line per GFA to phase, as
  `<gfa_filename> <start_edge>`. Here: `sample.gfa homA`.
- `mappings.txt` — synthetic barcode-to-contig k-mer mapping records in
  the `<read>_<barcode> <kmers> <contig>` format `main.cpp`'s
  `load_mappings` expects. Built to cover:
  - **BC001–BC008** (8 barcodes): map cleanly to `het1a` + `het2a` —
    should be called as one haplotype.
  - **BC009–BC013** (5 barcodes): map cleanly to `het1b` + `het2b` —
    the complementary haplotype.
  - **BC014, BC015**: map to only *one* bubble edge each — should be
    silently excluded (a barcode needs mappings to >1 edge to be
    usable, per `decide_barcode_haplotype_support`).
  - **BC016, BC017**: map to both bubbles but with a max per-edge score
    of 1 — should be excluded by the "noisy barcode" heuristic (max
    edge score must be > 1).
  - **BC001–BC004** also get a few k-mers against `homA`/`homB`
    (the homozygous backbone), to exercise `barcode_hom_mappings`.

## Running it

```
mkdir -p out
./phaser graphs.txt mappings.txt out
```

## Expected result

The run should report `1 graphs phased confidently`, with the winning
pair being `het1a,het2a` / `het1b,het2b`. `out/<prefix>.txt` should list
BC001–BC008 as supporting one haplotype and BC009–BC013 as supporting
the other, with BC014–BC017 absent from all three "barcodes supporting"
sections (they were filtered out before scoring).

Note: `main.cpp` currently derives the output filename by trimming the
GFA filename with `substr(start + 1, ...)`, which drops the first
character of the basename (e.g. `sample.gfa` → `out/ample.*`) — a
pre-existing bug, not specific to this dataset.
