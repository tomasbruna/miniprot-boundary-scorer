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
#define DEFAULT_EXON_SCORE 25
#define DEFAULT_INITIAL_EXON_SCORE 25
#define DEFAULT_INITIAL_INTRON_SCORE 0

void printUsage(char * name) {
    cout << "Usage: " << name << " < input -o output_file -s matrix_file "
            "[-w integer] [-k kernel] [-e min_exon_score] [-x min_initial_exon_score] [-i min_initial_intron_score] [-r]" << endl << endl;
    cout << "The program can parse multiple separate alignments saved in the same\n"
            "input. The input is read from stdin. Each input alignment is assumed\n"
            "to be on a single line (number of characters per line, controlled\n"
            "by -l option in Spaln, is larger than the alignment length)." << endl << endl;
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
            "      See the options -x and -i for details. Default = " <<
            DEFAULT_EXON_SCORE << endl;
    cout << "   -x Minimum initial exon score. Initial exons with lower scores\n"
            "      (as well as introns bordering such low-scoring exons and starts\n"
            "      inside them) are not printed. Initial exons with scores between\n"
            "      (-e and -x) must also define an initial intron which passes the\n"
            "      -i filter. Default = " <<
            DEFAULT_INITIAL_EXON_SCORE << endl;
    cout << "   -i Minimum initial intron score. Initial introns bordering\n"
            "      initial exons with scores < -e that have lower intron scores\n"
            "      (as well as initial exons bordering such low-scoring\n"
            "      introns and starts in those exons) are not printed.\n"
            "      Default = " <<
            DEFAULT_INITIAL_INTRON_SCORE << endl;
    cout << "   -r Process alignments on the reverse DNA strand (which are\n"
            "      ignored by default). This option might not be working properly\n"
            "      in this version!" << endl;
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
    bool processReverse = false;

    while ((opt = getopt(argc, argv, "o:w:s:k:e:i:x:r")) != EOF) {
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
            case 'r':
                processReverse = true;
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

    if (kernelType != "triangular" && kernelType != "box" &&
            kernelType != "parabolic" && kernelType != "triweight") {
        cerr << "error: Invalid kernel. Valid options are \"box\","
                "\"triangular\", \"parabolic\" and \"triweight\" kernels." << endl;
        printUsage(argv[0]);
        return 1;
    }

    ScoreMatrix * scoreMatrix = new ScoreMatrix();
    if (!scoreMatrix->loadFromFile(matrixFile)) {
        cerr << "error: Could not load scoring matrix" << endl;
        printUsage(argv[0]);
        return 1;
    }

    Kernel * kernel;
    if (kernelType == "triangular") {
        kernel = new TriangularKernel();
    } else if (kernelType == "box") {
        kernel = new BoxKernel();
    } else if (kernelType == "parabolic") {
        kernel = new ParabolicKernel();
    } else if (kernelType == "triweight") {
        kernel = new TriweightKernel();
    }

    Parser fileParser;
    fileParser.setWindowLegth(windowWidth);
    fileParser.setScoringMatrix(scoreMatrix);
    fileParser.setKernel(kernel);
    fileParser.setMinExonScore(minExonScore);
    fileParser.setMinInitialExonScore(minInitialExonScore);
    fileParser.setMinInitialIntronScore(minInitialIntronScore);
    fileParser.setProcessReverse(processReverse);

    int result = fileParser.parse(output);

    delete scoreMatrix;
    delete kernel;
    return result;
}

