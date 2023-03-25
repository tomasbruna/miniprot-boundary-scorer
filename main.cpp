#include "Parser.h"
#include "ScoreMatrix.h"
#include "Kernel.h"

#include <iostream>
#include <string>
#include <unistd.h>
#include <stdlib.h>

using namespace std;

#define DEFAULT_WINDOW_WIDTH 10
#define DEFAULT_KERNEL "triangular"
#define DEFAULT_EXON_SCORE -99999999999
#define DEFAULT_INITIAL_EXON_SCORE -99999999999
#define DEFAULT_INITIAL_INTRON_SCORE -99999999999
#define DEFAULT_GAP_PENALTY 4
#define DEFAULT_BOUNDARY_FRAMESHIFT_PENALTY 4
#define DEFAULT_EXON_FRAMEHSIFT_STOP_PENALTY 20

void printUsage(char * name) {
    cout << "Usage: " << name << " < input -o output_file -s matrix_file "
            "[-w integer] [-k kernel] [-e min_exon_score] [-x min_initial_exon_score] [-i min_initial_intron_score] [-g gap_penalty] [-f boundary_frameshift_penalty] [-F exon_frameshift_stop_penalty]" << endl << endl;
    cout << "The program parses the result of miniprot's \"--aln\" output. "
            "The input is read from stdin.\n" << endl << endl;
    cout << "Options:" << endl;
    cout << "   -o Where to save output file" << endl;
    cout << "   -s Path to amino acid scoring matrix" << endl;
    cout << "   -w Width of a scoring window around introns. Default = " <<
            DEFAULT_WINDOW_WIDTH << endl;
    cout << "   -k Specify type of weighting kernel used. Available opti-\n"
            "      ons are \"triangular\", \"box\", \"parabolic\" and \n"
            "      \"triweight\". Triangular kernel is the default option." << endl;
    cout << "   -e Minimum exon score. Exons with lower scores (as wells as int-\n"
            "      rons bordering such low-scoring exons and starts/stops inside\n"
            "      them are not printed). Initial exons are treated separately.\n"
            "      See the options -x and -i for details. Default = report all \n"
            "      exons." << endl;
    cout << "   -x Minimum initial exon score. Initial exons with lower scores\n"
            "      (as well as introns bordering such low-scoring exons and starts\n"
            "      inside them) are not printed. Initial exons with scores between\n"
            "      (-e and -x) must also define an initial intron which passes the\n"
            "      -i filter. Default = report all initial exons." << endl;
    cout << "   -i Minimum initial intron score. Initial introns bordering\n"
            "      initial exons with scores < -e that have lower intron scores\n"
            "      (as well as initial exons bordering such low-scoring\n"
            "      introns and starts in those exons) are not printed.\n"
            "      Default = report all initial introns." << endl;
    cout << "   -g Penalty for gaps, both in exons and around intron boundaries.\n"
            "      Default = " <<
            DEFAULT_GAP_PENALTY << endl;
    cout << "   -f Penalty for frameshifts around exon boundaries. After\n"
            "      a frameshift is detected, the rest of the exon boundary\n"
            "      is scored (still using the weighted score) by this penlalty,\n"
            "      regardless of the actual matches in the alignment.\n"
            "      Default = " <<
            DEFAULT_BOUNDARY_FRAMESHIFT_PENALTY << endl;
    cout << "   -F Penalty for frameshifts and read-through stop codons in \n"
            "      exons. Default = " <<
            DEFAULT_EXON_FRAMEHSIFT_STOP_PENALTY << endl;
}


int main(int argc, char** argv) {
    int opt;
    int windowWidth = DEFAULT_WINDOW_WIDTH;
    string output;
    string matrixFile = "";
    string kernelType = DEFAULT_KERNEL;
    double minExonScore = DEFAULT_EXON_SCORE;
    double minInitialIntronScore = DEFAULT_INITIAL_INTRON_SCORE;
    double minInitialExonScore = DEFAULT_INITIAL_EXON_SCORE;
    double gapPenalty = DEFAULT_GAP_PENALTY;
    double boundaryFrameshiftPenalty = DEFAULT_BOUNDARY_FRAMESHIFT_PENALTY;
    double exonFrameshiftStopPenalty = DEFAULT_EXON_FRAMEHSIFT_STOP_PENALTY;

    while ((opt = getopt(argc, argv, "o:w:s:k:e:i:x:f:g:F:")) != EOF) {
        switch (opt) {
            case 'o':
                output = optarg;
                break;
            case 'w':
                windowWidth = atoi(optarg);
                break;
            case 's':
                matrixFile = optarg;
                break;
            case 'k':
                kernelType = optarg;
                break;
            case 'e':
                minExonScore = atof(optarg);
                break;
            case 'i':
                minInitialIntronScore = atof(optarg);
                break;
            case 'x':
                minInitialExonScore = atof(optarg);
                break;
            case 'f':
                boundaryFrameshiftPenalty = atof(optarg);
                break;
            case 'F':
                exonFrameshiftStopPenalty = atof(optarg);
                break;
            case 'g':
                gapPenalty = atof(optarg);
                break;
            case '?':
                printUsage(argv[0]);
                return 1;
            default:
                return 1;
        }

    }

    if (output.size() == 0) {
        cerr << "error: Output file not specified" << endl;
        printUsage(argv[0]);
        return 1;
    }

    if (matrixFile.empty()) {
        cerr << "error: Score matrix not specified" << endl;
        printUsage(argv[0]);
        return 1;
    }

    if (minInitialExonScore > minExonScore) {
        cerr << "error: Minimum initial exon score must be lower than "
                "minimum exon score."<< endl;
        printUsage(argv[0]);
        return 1;
    }

    if (boundaryFrameshiftPenalty <= 0 || gapPenalty <= 0 ||
        exonFrameshiftStopPenalty <= 0) {
        cerr << "warning: The penalties are subtracted from the overall score "
                "and are thus assumed to be positive numbers. Perhaps "
                "specifying a non-positive penalty was an error?"<< endl;
    }

    ScoreMatrix * scoreMatrix = new ScoreMatrix();
    if (!scoreMatrix->loadFromFile(matrixFile)) {
        cerr << "error: Could not load scoring matrix" << endl;
        printUsage(argv[0]);
        return 1;
    }
    scoreMatrix->setPenalties(gapPenalty, boundaryFrameshiftPenalty,
                              exonFrameshiftStopPenalty);

    Kernel * kernel;
    if (kernelType == "triangular") {
        kernel = new TriangularKernel();
    } else if (kernelType == "box") {
        kernel = new BoxKernel();
    } else if (kernelType == "parabolic") {
        kernel = new ParabolicKernel();
    } else if (kernelType == "triweight") {
        kernel = new TriweightKernel();
    } else {
        cerr << "error: Invalid kernel. Valid options are \"box\","
                "\"triangular\", \"parabolic\" and \"triweight\" kernels." <<
                endl;
        printUsage(argv[0]);
        return 1;
    }

    Parser fileParser;
    fileParser.setWindowLegth(windowWidth);
    fileParser.setScoringMatrix(scoreMatrix);
    fileParser.setKernel(kernel);
    fileParser.setMinExonScore(minExonScore);
    fileParser.setMinInitialExonScore(minInitialExonScore);
    fileParser.setMinInitialIntronScore(minInitialIntronScore);

    int result = fileParser.parse(output);

    delete scoreMatrix;
    delete kernel;
    return result;
}

