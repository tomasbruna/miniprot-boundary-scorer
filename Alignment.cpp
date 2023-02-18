#include "Alignment.h"
#include "Parser.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <cmath>
#include <ctype.h>
#include <algorithm>

using namespace std;

Alignment::Alignment() {
    pairs.reserve(N);
    start = NULL;
    stop = NULL;
}

Alignment::~Alignment() {
    clear();
}

void Alignment::clear() {
    index = 0;
    insideIntron = false;
    donorFlag = false;
    introns.clear();
    for (unsigned int i = 0; i < exons.size(); i++) {
        delete exons[i];
    }
    exons.clear();
    dnaStart = 0;
    if (start != NULL) {
        delete start;
    }
    start = NULL;
    if (stop != NULL) {
        delete stop;
    }
    stop = NULL;
}

int Alignment::parse(istream& inputStream, string headerLine) {
    clear();
    vector<string> blockLines(BLOCK_ITEMS_CNT);
    string line;

    int status = parseHeader(headerLine);
    if (status != READ_SUCCESS) {
        cerr << "error: Invalid alignment header " << endl;
        return status;
    }

    // Load the alignment block and check its headers
    bool blockHeaderErr = false;
    for (i = 0; i < BLOCK_ITEMS_CNT; i++) {
        if (getline(inputStream, line) && !line.empty() && line.size() > BLOCK_OFFSET) {
            blockLines[i] = line;
            string lineStart = line.substr(0, BLOCK_OFFSET);
            switch(i) {
                case 0:
                    blockHeaderErr = (lineStart != "##ATN\t");
                    break;
                case 1:
                    blockHeaderErr = (lineStart != "##ATA\t");
                    break;
                case 2:
                    blockHeaderErr = (lineStart != "##AAS\t");
                    break;
                case 3:
                    blockHeaderErr = (lineStart != "##AQA\t");
                    break;
            }
            blockLines[i] = blockLines[i].substr(BLOCK_OFFSET);
        } else {
            blockHeaderErr = true;
        }
        if (blockHeaderErr) {
            cerr << "error: Unexpected alignment block headers" << endl;
            return FORMAT_FAIL;
        }
    }

    blockLength = blockLines[0].size();

    // Start actual parsing
    parseBlock(blockLines);

    // Potential single-exon gene
    if (exons.size() == 1) {
        exons.back()->initial = false;
    }

    return READ_SUCCESS;
}


int Alignment::parseHeader(string headerLine) {
    vector<string> cols;
    istringstream iss(headerLine);
    string col;
    while(getline(iss, col, '\t')) {
         cols.push_back(col);
    }
    if (cols.size() != 19) {
        cerr << "error: Unexpected number of columns in the header " << endl;
        return FORMAT_FAIL;
    }

    protein = cols[0];
    forward = false;
    if (cols[4] == "+") {
        forward = true;
    }
    // TODO: Load other relevant info.
    dnaStart = atoi(cols[7].c_str());
    realPositionCounter = dnaStart;
    return READ_SUCCESS;
}

void Alignment::parseBlock(const vector<string>& lines) {
    // Parse individual pairs
    for (unsigned i = 0; i < lines[0].size(); i++) {
        AlignedPair pair(lines[0][i], lines[1][i], lines[3][i]);

        checkForIntron(pair);
        checkForStart(pair);
        checkForStop(pair);

        if (pair.nucleotide != '-' ) {
            if (forward) {
                pair.realPosition = realPositionCounter++;
            } else {
                pair.realPosition = realPositionCounter--;
            }
        } else {
            // Don't increment the counter for gaps
            if (forward) {
                pair.realPosition = realPositionCounter - 1;
            } else {
                pair.realPosition = realPositionCounter + 1;
            }
        }

        // Reuse space if possible, just memory mgmt
        if ((int) pairs.size() <= index) {
            pairs.push_back(pair);
        } else {
            pairs[index] = pair;
        }
        index++;
    }
    assignCodonPhases();
}

bool Alignment::gapOrAA(char a) {
    if ((a >= 'A' && a <= 'Z') || a == '-') {
        return true;
    }
    return false;
}

bool Alignment::gapOrNT(char a) {
    if (a == 'A' || a == 'C' || a == 'G' || a == 'T' || a == 'N' ||
        a == 'a' || a == 'c' || a == 'g' || a == 't' || a == 'n' || a == '-') {
        return true;
    }
    return false;
}

void Alignment::assignCodonPhases() {
    int state = 0;
    for (unsigned int i = 0; i < blockLength; i++) {
        //TODO: What does a "$" sign mean? Probably a frameshift deletion.
        // Check for any other warnings of this type.
        if (gapOrAA(pairs[i].translatedCodon) != gapOrAA(pairs[i].protein)) {
            if ((pairs[i].translatedCodon != '*' || pairs[i].protein != '-') &&
                 i != blockLength - 3) {
                cerr << "Warning: Mismatch at alignment position " << i + 1 <<
                    " in the alignment of " << protein << endl;
            }
        }

        if (pairs[i].type != 'e') {
            continue;
        }

        if (gapOrAA(pairs[i].translatedCodon)) {
            state = 1;
        } else {
            state++;

            // Deal with out-of-phase alignment starts
            // TODO: Check whether this even happens.
            if (i == 0) {
                cerr << "Warning: and an out-of-phase alignment start" <<
                        " in the alignment of " << protein << endl;
                if (gapOrAA(pairs[i + 1].translatedCodon)) {
                    state = 3;
                } else {
                    state = 2;
                }
            }

            if (pairs[i].translatedCodon == '.'){
                pairs[i].translatedCodon = to_string(state).c_str()[0];
                pairs[i].protein = to_string(state).c_str()[0];
            }
        }
    }
}

void Alignment::checkForIntron(AlignedPair& pair) {
    int alignmentPosition;
    if (forward) {
        alignmentPosition = realPositionCounter - dnaStart + 1;
    } else {
        alignmentPosition = dnaStart - realPositionCounter + 1;
    }

    // Deal with the alignemnt start and end - a special case of intron boundary
    if (alignmentPosition == 1 || index == (int) blockLength - 1) {
        if (pair.type != 'e' || pair.nucleotide == '-') {
            cerr << "Warning: Unexpected alignment start/end" <<
                    " in the alignment of " << protein << endl;
        }

    }
    if (alignmentPosition == 1) {
        exons.push_back(new Exon(index));
    }
    if (index == (int) blockLength - 1) {
        exons.back()->end = index - 3;
    }

    if (donorFlag) {
        introns.back().donor[1] = pair.nucleotide;
        donorFlag = false;
    }

    if (!insideIntron && pair.type == 'i') { // intron start
        Intron i;
        i.start = index;
        i.donor[0] = pair.nucleotide;
        if (index != 0) {
            exons.back()->end = index - 1;
            i.leftExon = exons.back();
        }
        introns.push_back(i);
        insideIntron = true;
        donorFlag = true;
    } else if (insideIntron && pair.type != 'i') { // intron end
        insideIntron = false;

        if (pair.nucleotide != '-') {
            exons.push_back(new Exon(index));
        } else {
            exons.push_back(new Exon(index, true));
        }

        // Make the decision about exon phase based on how the preceeding exon
        // was split. Many checks remain from the Spaln parser that contained
        // strange exons, but it's OK to keep them as formatting checks.
        if (introns.back().start != 0) {
            if (gapOrAA(pairs[introns.back().start - 1].protein)) {
                exons.back()->phase = 2;
            } else if (introns.back().start - 2 < 0) {
                cerr << "Warning: Unexpected initial exon of length 1" <<
                        " in the alignment of " << protein << endl;
                exons.back()->phase = 1;
            } else {
                int position = introns.back().start - 2;
                if (pairs[position].type != 'e') {
                    // Check if still in exon (if the exon is just 1 nt long).
                    // If not, take last nt from next exon upstream
                    cerr << "Warning: Unexpected exon of length 1" <<
                        " in the alignment of " << protein << endl;
                    position = introns[introns.size() - 2].start - 1;
                }
                if (position >= 0 && gapOrAA(pairs[position].protein)) {
                    exons.back()->phase = 1;
                } else {
                    exons.back()->phase = 0;
                }
            }
        }

        introns.back().end = index - 1;
        introns.back().acceptor[0] = pairs[index - 2].nucleotide;
        introns.back().acceptor[1] = pairs[index - 1].nucleotide;
        // TODO: Another unneccessary check to remove
        if (introns.back().start != 0 && introns.back().gap == false) {
            introns.back().complete = true;
        }
        introns.back().rightExon = exons.back();
    }
}

void Alignment::checkForStart(AlignedPair& pair) {
    if (index == 2) {
        string codon = "";
        codon += pair.nucleotide;
        for (int i = 1; i < 3; i++) {
            codon = pairs[index - i].nucleotide + codon;
        }
        if (codon == "ATG") {
            // Check if protein alignment starts with its first M
            if (proteinStart == 1 && pairs[index - 1].protein == 'M') {
                start = new Codon(index - 2, exons.back());
                exons.back()->initial = true;
            }
        }
    }
}

void Alignment::checkForStop(AlignedPair& pair) {
    if (index == (int) blockLength - 1 && pairs[index - 3].type == 'e' &&
            (pair.nucleotide == 'a' || pair.nucleotide == 'g')) {
        string codon = "";
        codon += pair.nucleotide;
        for (int i = 1; i < 3; i++) {
            codon = pairs[index - i].nucleotide + codon;
        }

        if (codon == "taa" || codon == "tag" || codon == "tga") {
            stop = new Codon(index - 2, exons.back());
        }
    }
}

void Alignment::scoreHints(int windowWidth,
        const ScoreMatrix * scoreMatrix, Kernel * kernel) {
    this->scoreMatrix = scoreMatrix;
    this->kernel = kernel;
    this->kernel->setWidth(windowWidth);

    for (unsigned int i = 0; i < exons.size(); i++) {
        scoreExon(exons[i]);
    }

    scoreStart(windowWidth);
    scoreStop(windowWidth);

    for (unsigned int i = 0; i < introns.size(); i++) {
        if (introns[i].complete && !introns[i].scoreSet) {
            scoreIntron(introns[i], windowWidth);
        }
    }

}

double Alignment::scoreIntron(Intron& intron, int windowWidth) {
    intron.leftScore = intron.rightScore = 0;
    int left, right;

    // Determine if codon is split and how
    if (pairs[intron.start - 1].protein == '3' ||
            pairs[intron.start - 1].translatedCodon == '3') {
        // Codon is not split
        left = intron.start - 2;
        right = intron.end + 2;
    } else if (pairs[intron.start - 1].protein == '1'
            || pairs[intron.start - 1].translatedCodon == '1') {
        // Codon is split after the first nucleotide
        left = intron.start - 3;
        right = intron.end + 1;
    } else {
        // Codon is split after the second nucleotide
        left = intron.start - 1;
        right = intron.end + 3;
    }

    scoreLeft(intron, left, windowWidth);
    scoreRight(intron, right, windowWidth);
    double weightSum = kernel->weightSum();

    // Normalize alignments by the area under kernel
    if (intron.leftScore <= 0 || intron.rightScore <= 0) {
        intron.score = 0;
    } else {
        intron.score = (intron.leftScore / weightSum) *
                (intron.rightScore / weightSum);
        intron.score = sqrt(intron.score);
    }

    intron.score /= scoreMatrix->getMaxScore();
    intron.scoreSet = true;
    return intron.score;
}

void Alignment::scoreLeft(Intron & intron, int start, int windowWidth) {
    for (int i = start; i > (start - windowWidth * 3); i -= 3) {
        // Check for end of local alignment
        if (i < 0 || pairs[i].type != 'e') {
            return;
        }
        double weight = kernel->getWeight((i - start) / 3);
        intron.leftScore += pairs[i].score(scoreMatrix) * weight;
    }
}

void Alignment::scoreRight(Intron & intron, int start, int windowWidth) {
    for (int i = start; i < (start + windowWidth * 3); i += 3) {
        // Check for end of local alignment
        if (i >= index || pairs[i].type != 'e') {
            return;
        }
        double weight = kernel->getWeight((i - start) / 3);
        intron.rightScore += pairs[i].score(scoreMatrix) * weight;
    }
}

void Alignment::scoreStart(int windowWidth) {
    if (start == NULL) {
        return;
    }
    start->score = 0;
    for (int i = start->position + 1; i < (start->position + 1 + windowWidth * 3); i += 3) {
        // Check for end of local alignment
        if (i >= index || pairs[i].type != 'e') {
            break;
        }
        double weight = kernel->getWeight((i - start->position - 1) / 3);
        start->score += pairs[i].score(scoreMatrix) * weight;
    }
    start->score /=  kernel->weightSum();
    start->score /= scoreMatrix->getMaxScore();
}

void Alignment::scoreStop(int windowWidth) {
    if (stop == NULL) {
        return;
    }
    stop->score = 0;
    for (int i = stop->position - 2; i > (stop->position - 2 - windowWidth * 3); i -= 3) {
        // Check for end of local alignment
        if (i < 0 || pairs[i].type != 'e') {
            break;
        }
        double weight = kernel->getWeight((i - stop->position + 2) / 3);
        stop->score += pairs[i].score(scoreMatrix) * weight;
    }
    stop->score /=  kernel->weightSum();
    stop->score /= scoreMatrix->getMaxScore();
}

void Alignment::scoreExon(Exon* exon) {
    int i = exon->start;
    exon->score = 0;
    int length = 0;
    while (i <= exon->end) {
        if (gapOrAA(pairs[i].protein)) {
            exon->score += pairs[i].score(scoreMatrix);
            length++;
        }
        i++;
    }
    exon->normalizedScore = exon->score / length;
}

void Alignment::printHints(string output, double minExonScore,
                           double minInitialExonScore, double minInitialIntronScore) {
    ofstream ofs(output.c_str(), std::ofstream::out | std::ofstream::app);
    char strand;
    if (forward) {
        strand = '+';
    } else {
        strand = '-';
    }

    printIntrons(ofs, strand, minExonScore, minInitialExonScore,
                 minInitialIntronScore);
    printStart(ofs, strand, minExonScore, minInitialExonScore,
               minInitialIntronScore);
    printExons(ofs, strand, minExonScore, minInitialExonScore,
               minInitialIntronScore);
    printStop(ofs, strand, minExonScore);

    ofs.close();
}

void Alignment::printIntrons(ofstream& ofs, char strand,
                             double minExonScore, double minInitialExonScore,
                             double minInitialIntronScore) {
    for (unsigned int i = 0; i < introns.size(); i++) {
        if (!introns[i].complete || introns[i].rightExon->score < minExonScore) {
            continue;
        }

        if (introns[i].leftExon->score < minExonScore) {
            // Initial intron can have lower leftExon score, if it passes the
            // intron score filter
            if (!introns[i].leftExon->initial ||
                introns[i].leftExon->score < minInitialExonScore ||
                introns[i].score < minInitialIntronScore) {
                continue;
            }
        }

        string spliceSites(introns[i].donor, 2);
        spliceSites.append("_");
        spliceSites.append(introns[i].acceptor, 2);

        ofs << gene << "\tminiprot_scorer\tIntron\t";
        if (forward) {
            ofs << pairs[introns[i].start].realPosition << "\t";
            ofs << pairs[introns[i].end].realPosition << "\t";
        } else {
            ofs << pairs[introns[i].end].realPosition << "\t";
            ofs << pairs[introns[i].start].realPosition << "\t";
        }
        ofs << ".\t" << strand << "\t.\tprot=" << protein;
        ofs << "; intron_id=" << i + 1 << ";";
        ofs << " initial=" << introns[i].leftExon->initial << ";";
        ofs << " splice_sites=" << spliceSites << ";";
        ofs << " al_score=" << introns[i].score << ";";
        ofs << " LeScore=" << introns[i].leftExon->score << ";";
        ofs << " ReScore=" << introns[i].rightExon->score << ";";
        ofs << " LeNScore=" << introns[i].leftExon->normalizedScore << ";\n";
    }
}

void Alignment::printStart(ofstream& ofs, char strand,
                           double minExonScore,
                           double minInitialExonScore,
                           double minInitialIntronScore) {
    if (start == NULL || start->exon->score < minInitialExonScore) {
        return;
    }

    // Print start with low exon score if the next intron passes filters
    if (start->exon->score < minExonScore) {
        if (introns.size() == 0 || !introns[0].complete ||
            !introns[0].leftExon->initial ||
            introns[0].rightExon->score < minExonScore ||
            introns[0].score < minInitialIntronScore) {
            return;
        }
    }

    ofs << gene << "\tminiprot_scorer\tstart_codon\t";
    if (forward) {
        ofs << pairs[start->position].realPosition << "\t";
        ofs << pairs[start->position + 2].realPosition  << "\t";
    } else {
        ofs << pairs[start->position + 2].realPosition  << "\t";
        ofs << pairs[start->position].realPosition << "\t";
    }
    ofs << ".\t" << strand << "\t0\tprot=" << protein << ";";
    ofs << " al_score=" << start->score << ";";
    ofs << " eScore=" << start->exon->score << ";";
    ofs << " eNScore=" << start->exon->normalizedScore << ";";

    // Only save next intron coordinates if the intron passes filters
    if (introns.size() != 0 && introns[0].complete &&
        introns[0].rightExon->score >= minExonScore &&
        introns[0].score >= minInitialIntronScore &&
        introns[0].leftExon->initial) {
        int offsetStart = pairs[introns[0].start].realPosition -
            pairs[start->position].realPosition;
        int offsetEnd = pairs[introns[0].end].realPosition -
            pairs[start->position].realPosition;
        ofs << " nextIntron=" << offsetStart << "-" << offsetEnd << ";\n";
    } else {
        ofs << " nextIntron=-;\n";
    }
}

void Alignment::printExons(ofstream& ofs, char strand, double minExonScore,
                           double minInitialExonScore,
                           double minInitialIntronScore) {
    for (unsigned int i = 0; i < exons.size(); i++) {
        if (exons[i]->score < minExonScore) {
            // Initial exons can have lower score, if the first intron
            // passes filters
            if (exons[i]->score < minInitialExonScore || !exons[i]->initial ||
                start->score <= 0 || introns.size() == 0 ||
                !introns[0].leftExon->initial || !introns[0].complete ||
                introns[0].rightExon->score < minExonScore ||
                introns[0].score < minInitialIntronScore) {
                continue;
            }
        }
        ofs << gene << "\tminiprot_scorer\tCDS\t";
        if (forward) {
            if (exons[i]->gapStart) {
                ofs << pairs[exons[i]->start].realPosition + 1 << "\t";
            } else {
                ofs << pairs[exons[i]->start].realPosition << "\t";
            }
            ofs << pairs[exons[i]->end].realPosition << "\t";
        } else {
            ofs << pairs[exons[i]->end].realPosition << "\t";
            ofs << pairs[exons[i]->start].realPosition << "\t";
        }
        ofs << ".\t" << strand << "\t" << exons[i]->phase << "\tprot=" << protein;
        ofs << "; exon_id=" << i + 1 << ";";
        ofs << " initial=" << exons[i]->initial << ";";
        ofs << " eScore=" << exons[i]->score << ";";
        ofs << " eNScore=" << exons[i]->normalizedScore << ";\n";
    }
}

void Alignment::printStop(ofstream& ofs, char strand, double minExonScore) {
    if (stop != NULL && stop->exon->score >= minExonScore) {
        ofs << gene << "\tminiprot_scorer\tstop_codon\t";
        if (forward) {
            ofs << pairs[stop->position].realPosition << "\t";
            ofs << pairs[stop->position + 2].realPosition  << "\t";
        } else {
            ofs << pairs[stop->position + 2].realPosition  << "\t";
            ofs << pairs[stop->position].realPosition << "\t";
        }
        ofs << ".\t" << strand << "\t0\tprot=" << protein << ";";
        ofs << " al_score=" << stop->score << ";";
        ofs << " eScore=" << stop->exon->score << ";\n";
    }
}

void Alignment::print(ostream& os) {
    for (unsigned int i = 0; i < blockLength; i++) {
        os << pairs[i].nucleotide;
    }
    os << endl;
    for (unsigned int i = 0; i < blockLength; i++) {
        os << pairs[i].translatedCodon;
    }
    os << endl;
    for (unsigned int i = 0; i < blockLength; i++) {
        os << pairs[i].type;
    }
    os << endl;
    for (unsigned int i = 0; i < blockLength; i++) {
        os << pairs[i].protein;
    }
    os << endl;
}

string Alignment::getGene() {
    return gene;
}

string Alignment::getProtein() {
    return protein;
}

int Alignment::getLength() {
    return index;
}

Alignment::AlignedPair::AlignedPair(char n, char tc, char p) :
nucleotide(n),
translatedCodon(tc),
protein(p) {

    if (!gapOrNT(n)) {
        cerr << "Warning: Unexpected nucleotide \"" << n << "\". " <<
                "This is probably miniprot's compact notation for long "
                "introns. Not supported yet by this parser." << endl;
    }

    if (n == '-') {
        this->type = 'e';
    } else if (islower(n)) {
        this->type = 'i';
    } else {
        this->type = 'e';
    }

    // TODO: Dealing with the stop...
    if (translatedCodon == '*') {
        translatedCodon = 'A';
    }
}

double Alignment::AlignedPair::score(const ScoreMatrix * scoreMatrix) {
    return scoreMatrix->getScore(translatedCodon, protein);
}

Alignment::Intron::Intron() {
    scoreSet = false;
    complete = false;
    gap = false;
}

Alignment::Exon::Exon(int start, bool gapStart) {
    this->start = start;
    this->gapStart = gapStart;
    phase = 0;
    scoreSet = false;
    initial = false;
}

Alignment::Codon::Codon(int position, Exon * exon) {
    this->position = position;
    this->exon = exon;
}
