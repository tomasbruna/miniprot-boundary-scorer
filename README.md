# Miniprot boundary scorer

Miniprot boundary scorer parses introns, starts, stops and exons from [miniprot's](https://github.com/lh3/miniprot) alignment output and scores them. Introns, starts and stops are scored based on local alignment quality around their boundaries. Detailed description of how the scores are computed is available in https://academic.oup.com/nargab/article/2/2/lqaa026/5836691.

## Installation and usage

To install, clone this repository and run `make` in the root folder.

To run, use the following command:

    miniprot_boundary_scorer < miniprot_input -o output_file -s matrix_file [-w integer] [-k kernel] [-e min_exon_score]

Input details:

* The program can parse multiple separate alignments saved in the same input.
* The input is read from stdin.

Available options are:

```
   -o Where to save output file
   -s Path to amino acid scoring matrix
   -w Width of a scoring window around introns. Default = 10
   -k Specify type of weighting kernel used. Available opti-
      ons are "triangular", "box", "parabolic" and 
      "triweight". Triangular kernel is the default option.
   -e Minimum exon score. Exons with lower scores (and int-
      rons bordering such exons; start and stops inside the 
      exons) are not printed. Default = 25
```

## Tests

Unit tests are located in the `test` folder. To compile a test binary, run
`make test` in the root folder. Subsequently, run the test binary to evaluate
the tests:

    make test
    test/t_miniprot_boundary_scorer
