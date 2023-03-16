#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <fstream>
#include <vector>
#include <string>
#include "ScoreMatrix.h"
#include "Kernel.h"

using namespace std;

/// Class for parsing a single DNA-protein alignment

class Alignment {
public:
    Alignment();
    ~Alignment();
    /**
     * Parse a single DNA-protein alignment.
     * The function checks if the general structure of the alignment
     * is OK but it does not check the validity of every single base/protein.
     *
     * @param fstream    File stream starting at the position of
     *                    the alignment start.
     * @param headerLine Header line associated with this alignment
     */
    int parse(istream & inputStream, string headerLine);
    /**
     * @return Name of the aligned seqid
     */
    string getSeqid();
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
        AlignedPair(char n, char tc, char p);
        /**
         * Return amino acid score
         *
         * The intron frameshift flag indicates that the intron border is in
         * a "frameshift" state and should be penalized with the
         * INTRON_FRAMESHIFT penalty instead.
         *
         * @return AA score
         */
        double score(const ScoreMatrix * scoreMatrix,
                     bool introFrameshift = false);
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
        Exon(int start, bool gapStart = false);
        int start, end;
        double score;
        double normalizedScore;
        int phase;
        bool scoreSet;
        bool initial;
        bool gapStart;
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
        double leftScore, rightScore;
        double printedLength;
        double realLength;
        /// Flag indicating that a gap, or aligned protein, was detected
        /// inside intron
        Exon * leftExon;
        Exon * rightExon;
    };

    /**
     * Clear the object for a new alignment pair
     */
    void clear();
    /**
     *  Parse the alignment header
     */
    int parseHeader(string headerLine);
    /**
     *  Get a named attribute from the header
     */
    string getHeaderAttribute(const vector<string>& cols, string attribute);
    /**
     *  Check for the intron compact notation and parse it
     */
    bool checkForCompactNotation(char a);
    /**
     *  Parse individual block of lines containing the alignment and its properties
     */
    void parseBlock(const vector<string>& lines);
    /**
     * Check if the given character is an amino acid, stop codon, or gap
     */
    bool gapStopOrAA(char a);
    /**
     * Check if the given character represents frameshift
     */
    bool frameshift(char a);
    /**
     * Check if the given character is a nucleotide or a gap
     */
    static bool gapOrNT(char a);
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

    static const int BLOCK_ITEMS_CNT = 4;
    static const int BLOCK_OFFSET = 6;
    /// Starting position of the alignment in DNA
    int dnaStart;
    /// Starting position of the alignment in protein
    int proteinStart;
    int proteinLength;
    int alignedProteinLength;
    int AS;
    int ms;
    int positiveMatches;
    int exactMatches;
    int i;
    string seqid;
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
    bool  compactState;
    string compactLength;
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


