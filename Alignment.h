#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <fstream>
#include <vector>
#include <string>
#include "ScoreMatrix.h"
#include "Kernel.h"

using namespace std;

/// Class for parsing a single gene-protein alignment

class Alignment {
public:
    Alignment();
    ~Alignment();
    /**
     * Parse a single gene-protein alignment.
     * The function checks if the general structure of the alignment
     * is OK but it does not check the validity of every single base/protein.
     *
     * @param fstream    File stream starting at the position of
     *                    the alignment start.
     * @param headerLine Header line associated with this alignment
     * @param forward    Whether the alignment is on forward strand
     */
    int parse(istream & inputStream, string headerLine, bool forward);
    /**
     * @return Name of the aligned gene
     */
    string getGene();
    /**
     * @return Name of the aligned protein
     */
    string getProtein();
    /**
     * @return Total alignment length in nucleotides
     */
    int getLength();
    /**
     * Print the whole alignment
     */
    void print(ostream & os);
    /**
     * Print scored hints
     *
     * @param output                Output file name
     * @param minExonScore          Do not print hints with exon score lower than this
     * @param minInitialExonScore   Do not print hints with initial exon score lower than this
     * @param minInitialIntronScore Do not print hints with initial intron score lower than this
     */
    void printHints(string output, double minExonScore,
                    double minInitialExonScore, double minInitialIntronScore);
    /**
     * Score all hints in the alignment
     * @param windowWidth Number of amino acids scored in the upstream/
     *                    downstream regions
     * @param scoreMatrix Scoring matrix used for scoring amino acids.
     */
    void scoreHints(int windowWidth,
            const ScoreMatrix * scoreMatrix, Kernel * kernel);
private:
    /// Single nucleotide-amino acid pair
    struct AlignedPair {
        /**
         * Save pair and determine exon/intron
         */
        AlignedPair(char tc, char n, char p, bool insideIntron);
        /**
         * Return amino acid score
         * @return AA score
         */
        double score(const ScoreMatrix * scoreMatrix);
        char nucleotide;
        /**
         * The protein translations and proteins are saved as follows: 1A3
         * Where A is the AA, numbers 1 and 3 fill
         * the space created by 3 nucleotide to 1 AA mapping.
         */
        char translatedCodon;
        char protein;
        /**
         * Type of DNA base (based on alignment)
         * 'i' for intron
         * 'e' for exon
         */
        char type;
        /// Position of a nucleotide in the alignment relative to a gene start
        int realPosition;
    };

    /// Structure for parsed exons
    struct Exon {
        Exon(int start);
        int start, end;
        double score;
        double normalizedScore;
        int phase;
        bool scoreSet;
        bool initial;
    };

    /// Structure for parsed starts and stops
    struct Codon {
        Codon(int position, Exon * exon);
        int position;
        double score;
        Exon * exon;
    };

    /// Structure for parsed introns
    struct Intron {
        Intron();
        unsigned int start, end;
        double score;
        bool scoreSet;
        char donor[2];
        char acceptor[2];
        bool complete;
        double leftScore, rightScore;
        /// Flag indicating that a gap, or aligned protein, was detected
        /// inside intron
        bool gap;
        Exon * leftExon;
        Exon * rightExon;
    };

    /**
     * Clear the object for a new alignment pair
     */
    void clear();
    /**
     *  Parse gene and protein name from header line
     */
    int parseHeader(string headerLine);
    /**
     *  Parse individual block of lines containing the alignment and its properties
     */
    void parseBlock(const vector<string>& lines);
    /**
     * Check if the given character is an amino acid or a gap
     */
    bool gapOrAA(char a);
    /**
     * Assign phases to all positions
     */
    void assignCodonPhases();
    void printLineError();
    /**
     * Detect and save introns.
     * The function also retrieves information associated with the intron
     * such as its start/end position and donor/acceptor site.
     * Exons are saved as well during this process.
     */
    void checkForIntron(AlignedPair & pair);
    /**
     * Detect and save start codon
     */
    void checkForStart(AlignedPair & pair);
    /**
     * Detect and save stop codon
     */
    void checkForStop(AlignedPair & pair);
    /**
     * Determine score of a single intron using exon alignment in the
     * upstream and downstream region
     */
    double scoreIntron(Intron & intron, int windowWidth);
    /**
     * Compute alignment score of amino acids upstream of intron
     */
    void scoreLeft(Intron & intron, int start, int windowWidth);
    /**
     * Compute alignment score of amino acids downstream of intron
     */
    void scoreRight(Intron & intron, int start, int windowWidth);
    void scoreExon(Exon * exon);
    void scoreStart(int windowWidth);
    void scoreStop(int windowWidth);

    void printIntrons(ofstream & ofs, char strand, double minExonScore,
                      double minInitialExonScore, double minInitialIntronScore);
    void printStart(ofstream & ofs, char strand, double minExonScore,
                    double minInitialExonScore, double minInitialIntronScore);
    void printExons(ofstream & ofs, char strand, double minExonScore,
                    double minInitialExonScore, double minInitialIntronScore);
    void printStop(ofstream & ofs, char strand, double minExonScore);

    static const int BLOCK_ITEMS_CNT = 3;
    static const int BLOCK_OFFSET = 9;
    /// Starting position of the alignment in DNA
    int dnaStart;
    /// Starting position of the alignment in protein
    int proteinStart;
    int i;
    string gene;
    string protein;
    /// Overall alignment length
    unsigned int blockLength;
    /// Overall position in alignment, including gaps
    int index;
    /// Track position of nucleotides in the alignment relative to seed start
    /// (gaps do not increment the counter)
    int realPositionCounter;
    bool forward;
    /// Initial size of alignment vector
    static const int N = 3000;
    /// Array containing the actual alignment pairs
    vector<AlignedPair> pairs;
    // Whether the parser is inside intron state
    bool insideIntron;
    /// Flag indicating that donor position of an intron is being read
    bool donorFlag;
    vector<Intron> introns;
    vector<Exon*> exons;
    Codon * start;
    Codon * stop;
    const ScoreMatrix * scoreMatrix;
    Kernel * kernel;
};


#endif /* ALIGNMENT_H */


