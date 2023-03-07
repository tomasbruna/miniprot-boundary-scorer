# Miniprot boundary scorer

Miniprot boundary scorer parses introns, starts, stops and exons from [miniprot's](https://github.com/lh3/miniprot) alignment output and scores them. Introns, starts and stops are scored based on local alignment quality around their boundaries. The scoring is similar to [ProtHint's](https://github.com/gatech-genemark/ProtHint) scoring of Spaln alignments, described in https://academic.oup.com/nargab/article/2/2/lqaa026/5836691.

## Installation and usage

To install, clone this repository and run `make` in the root folder.

To run, use the following command:

    miniprot_boundary_scorer < miniprot_input -o output_file -s matrix_file [optional arguments]

Input details:

* The program can parse multiple separate alignments saved in the same input.
* The input is read from stdin -- can be piped directly from miniprot.

Available options are:

```
   -o Where to save output file
   -s Path to amino acid scoring matrix
   -w Width of a scoring window around introns. Default = 10
   -k Specify type of weighting kernel used. Available opti-
      ons are "triangular", "box", "parabolic" and 
      "triweight". Triangular kernel is the default option.
   -e Minimum exon score. Exons with lower scores (as wells as int-
      rons bordering such low-scoring exons and starts/stops inside
      them are not printed). Initial exons are treated separately.
      See the options -x and -i for details. Default = 25
   -x Minimum initial exon score. Initial exons with lower scores
      (as well as introns bordering such low-scoring exons and starts
      inside them) are not printed. Initial exons with scores between
      (-e and -x) must also define an initial intron which passes the
      -i filter. Default = 25
   -i Minimum initial intron score. Initial introns bordering
      initial exons with scores < -e that have lower intron scores
      (as well as initial exons bordering such low-scoring
      introns and starts in those exons) are not printed.
      Default = 0
   -g Penalty for gaps, both in exons and around intron boundaries.
      Default = 4
   -f Penalty for frameshifts around exon boundaries. After
      a frameshift is detected, the rest of the exon boundary
      is scored (still using the weighted score) by this penlalty,
      regardless of the actual matches in the alignment.
      Default = 4
   -F Penalty for frameshifts and read-through stop codons in
      exons. Default = 20
```

## Tests

Unit tests are located in the `test` folder. To compile a test binary, run
`make test` in the root folder. Subsequently, run the test binary to evaluate
the tests:

    make test
    test/t_miniprot_boundary_scorer
