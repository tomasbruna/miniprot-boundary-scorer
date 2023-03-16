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
    alignedProteinLength = 0;
    AS = 0;
    ms = 0;
    positiveMatches = 0;
    exactMatches = 0;
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
    compactState = false;
    compactLength = "";
}

int Alignment::parse(istream& inputStream, string headerLine) {
    clear();
    vector<string> blockLines(BLOCK_ITEMS_CNT);
    string line;

    int status = parseHeader(headerLine);
    if (status != READ_SUCCESS) {
        if (status != EMPTY_ALIGNMENT) {
            cerr << "error: Invalid alignment header " << endl;
        }
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
            blockLines[i] += "<";
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

string Alignment::getHeaderAttribute(const vector<string>& cols, string attribute) {
    for (i=0; i < cols.size(); i++) {
        if (cols[i].substr(0, attribute.size()) == attribute) {
            return cols[i].substr(attribute.size() + 1);
        }
    }
    cerr << "warning: attribute " << attribute << " was not found in the "
         << "header. Returning 0 instead" << endl;
    return "0";
}

int Alignment::parseHeader(string headerLine) {
    vector<string> cols;
    istringstream iss(headerLine);
    string col;
    while(getline(iss, col, '\t')) {
         cols.push_back(col);
    }
    if (cols.size() == 12) {
        return EMPTY_ALIGNMENT;
    }

    if (cols.size() < 19) {
        cerr << "error: Unexpected number of columns in the header. Miniprot "
             << "PAF header should contain 19 fields" << endl;
        return FORMAT_FAIL;
    }

    protein = cols[0];
    proteinLength = atoi(cols[1].c_str());
    proteinStart = atoi(cols[2].c_str()) + 1;

    forward = false;
    // miniprot is using bed-style coordinates
    if (cols[4] == "+") {
        forward = true;
        dnaStart = atoi(cols[7].c_str()) + 1;
    } else {
        dnaStart = atoi(cols[8].c_str());
    }
    realPositionCounter = dnaStart;
    seqid = cols[5].c_str();
    AS = atoi(getHeaderAttribute(cols, "AS:i").c_str());
    ms = atoi(getHeaderAttribute(cols, "ms:i").c_str());
    return READ_SUCCESS;
}

bool Alignment::checkForCompactNotation(char a) {
    if (a == '~') {
        if (!insideIntron) {
            cerr << "Warning: Unexpected character \"~\" outside of " <<
                    "an intron in the alignment of " << protein << endl;
        }
        if (compactState) {
            if (introns.back().realLength != -1) {
                cerr << "Warning: Unexpected compact intron length " <<
                        "representation in the alignment of " << protein << endl;
            }
            introns.back().realLength = atoi(compactLength.c_str());
            compactLength = "";
            compactState = false;
        } else {
            compactState = true;
        }
        return true;
    }
    if (compactState) {
        if (a < '0' || a > '9') {
            cerr << "Warning: A number expected after \"~\" " <<
                    "in the alignment of " << protein << endl;
        }
        compactLength += a;
        return true;
    }
    return false;
}

void Alignment::parseBlock(const vector<string>& lines) {
    // Parse individual pairs
    for (unsigned i = 0; i < lines[0].size(); i++) {

        if (checkForCompactNotation(lines[0][i])) {
            continue;
        }

        AlignedPair pair(lines[0][i], lines[1][i], lines[3][i]);

        // Format checks
        if (gapStopOrAA(pair.translatedCodon) != gapStopOrAA(pair.protein)) {
            if (i + 3 < blockLength && lines[1][i + 3] == '<') {
                // Alignment end
            } else if (!frameshift(pair.translatedCodon)) {
                cerr << "Warning: Mismatch at alignment position " << i + 1 <<
                    " in the alignment of " << protein << endl;
            }
        }

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

        if (lines[3][i] >= 'A' && lines[3][i] <= 'Z') {
            alignedProteinLength++;
        }

        // Reuse space if possible, just memory mgmt
        if ((int) pairs.size() <= index) {
            pairs.push_back(pair);
        } else {
            pairs[index] = pair;
        }
        index++;
    }
}

bool Alignment::gapStopOrAA(char a) {
    if ((a >= 'A' && a <= 'Z') || a == '-' || a == '*') {
        return true;
    }
    return false;
}

bool Alignment::frameshift(char a) {
    if (a == '$' || a == '!') {
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

void Alignment::checkForIntron(AlignedPair& pair) {
    int alignmentPosition;
    if (forward) {
        alignmentPosition = realPositionCounter - dnaStart + 1;
    } else {
        alignmentPosition = dnaStart - realPositionCounter + 1;
    }

    // Deal with the alignemnt start and end - a special case of intron boundary
    if (alignmentPosition == 1 || pair.nucleotide == '<') {
        if (pair.type != 'e' || pair.nucleotide == '-') {
            cerr << "Warning: Unexpected alignment start/end" <<
                    " in the alignment of " << protein << endl;
        }

    }
    if (alignmentPosition == 1) {
        exons.push_back(new Exon(index));
    }
    if (pair.nucleotide == '<') {
        // Proteins ending with a * in the fasta input.
        // TODO: To properly test this, create a protein input where
        // every protein ends with a stop codon
        if (pairs[index - 6].protein == '*') {
            exons.back()->end = index - 7;
        } else {
            exons.back()->end = index - 4;
        }
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
        if (gapStopOrAA(pairs[introns.back().start - 1].protein)) {
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
            if (position >= 0 && gapStopOrAA(pairs[position].protein)) {
                exons.back()->phase = 1;
            } else {
                exons.back()->phase = 0;
            }
        }

        introns.back().end = index - 1;
        introns.back().acceptor[0] = pairs[index - 2].nucleotide;
        introns.back().acceptor[1] = pairs[index - 1].nucleotide;
        introns.back().rightExon = exons.back();

        // Add the length stored in the compact notation
        if (introns.back().realLength != -1) {
            int compactOffset = introns.back().realLength -
                                (introns.back().end - introns.back().start + 1);
            if (forward) {
                realPositionCounter += compactOffset;
                // This is kind of hacky, but only the boundary coordinate
                // needs to be correct.
                pairs[index - 1].realPosition = realPositionCounter - 1;
            } else {
               realPositionCounter -= compactOffset;
               pairs[index - 1].realPosition = realPositionCounter + 1;
            }
        }
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
            if (proteinStart == 1 && pairs[index - 2].protein == 'M') {
                start = new Codon(index - 2, exons.back());
                exons.back()->initial = true;
            }
        }
    }
}

void Alignment::checkForStop(AlignedPair& pair) {
    if (pair.nucleotide == '<') {
        if (pairs[index - 3].translatedCodon == '*') {
            stop = new Codon(index - 3, exons.back());
        }
    }
    if (pair.protein == '*') {
        cerr << "Warning: A stop codon detected inside protein " <<
                protein << " at al. position " << index << endl;
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
        if (!introns[i].scoreSet) {
            scoreIntron(introns[i], windowWidth);
        }
    }

}

double Alignment::scoreIntron(Intron& intron, int windowWidth) {
    intron.leftScore = intron.rightScore = 0;
    int left, right;

    // Determine if a codon is split and how
    // A separate check for upstream and downstream boundaries is needed
    // because of a possibility of rare frameshifts at the intron boundary

    // is Frameshift function
    char c1 = pairs[intron.start - 1].translatedCodon;
    char c2 = pairs[intron.start - 2].translatedCodon;
    if (gapStopOrAA(c1) || frameshift(c1)) {
        left = intron.start - 1;
    } else if (gapStopOrAA(c2) || frameshift(c2)) {
        left = intron.start - 2;
    } else {
        left = intron.start - 3;
    }

    c1 = pairs[intron.end + 1].translatedCodon;
    c2 = pairs[intron.end + 2].translatedCodon;
    if (gapStopOrAA(c1) || frameshift(c1)) {
        right = intron.end + 1;
    } else if (gapStopOrAA(c2) || frameshift(c2)) {
        right = intron.end + 2;
    } else {
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
    bool fs = false;
    for (int i = start; i > (start - windowWidth * 3); i -= 1) {
        // Check for end of local alignment
        if (i < 0 || pairs[i].type != 'e') {
            return;
        }
        if (frameshift(pairs[i].translatedCodon)) {
            fs = true;
        }
        if (i % 3 == start % 3) {
            double weight = kernel->getWeight((i - start) / 3);
            intron.leftScore += pairs[i].score(scoreMatrix, fs) * weight;
        }
    }
}

void Alignment::scoreRight(Intron & intron, int start, int windowWidth) {
    bool fs = false;
    for (int i = start; i < (start + windowWidth * 3); i += 1) {
        // Check for end of local alignment.
        if (i >= exons.back()->end + 1 || pairs[i].type != 'e') {
            return;
        }
        if (frameshift(pairs[i].translatedCodon)) {
            fs = true;
        }
        if (i % 3 == start % 3) {
            double weight = kernel->getWeight((i - start) / 3);
            intron.rightScore += pairs[i].score(scoreMatrix, fs) * weight;
        }
    }
}

void Alignment::scoreStart(int windowWidth) {
    if (start == NULL) {
        return;
    }
    start->score = 0;
    bool fs = false;
    for (int i = start->position; i < (start->position + windowWidth * 3); i += 1) {
        // Check for end of local alignment
        if (i >= exons.back()->end + 1  || pairs[i].type != 'e') {
            break;
        }
        if (frameshift(pairs[i].translatedCodon)) {
            fs = true;
        }
        if (i % 3 == start->position % 3) {
            double weight = kernel->getWeight((i - start->position) / 3);
            start->score += pairs[i].score(scoreMatrix, fs) * weight;
        }
    }
    start->score /=  kernel->weightSum();
    start->score /= scoreMatrix->getMaxScore();
}

void Alignment::scoreStop(int windowWidth) {
    if (stop == NULL) {
        return;
    }
    stop->score = 0;
    bool fs = false;
    for (int i = stop->position - 3; i > (stop->position - 3 - windowWidth * 3); i -= 1) {
        // Check for end of local alignment.
        if (i < 0 || pairs[i].type != 'e') {
            break;
        }
        if (frameshift(pairs[i].translatedCodon)) {
            fs = true;
        }
        if (i % 3 == (stop->position - 3) % 3) {
            double weight = kernel->getWeight((i - stop->position + 3) / 3);
            stop->score += pairs[i].score(scoreMatrix, fs) * weight;
        }
    }
    stop->score /=  kernel->weightSum();
    stop->score /= scoreMatrix->getMaxScore();
}

void Alignment::scoreExon(Exon* exon) {
    int i = exon->start;
    exon->score = 0;
    int length = 0;
    while (i <= exon->end) {
        if (pairs[i].translatedCodon == '!' || pairs[i].translatedCodon == '$') {
            // A "single" frameshift can be represented by more than one symbol
            // depending on its phase. Penalize it only once.
            if (pairs[i - 1].translatedCodon != '!' && pairs[i - 1].translatedCodon != '$') {
                exon->score += pairs[i].score(scoreMatrix);
                length++;
            }
        } else if (gapStopOrAA(pairs[i].translatedCodon)) {
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
        if (introns[i].rightExon->score < minExonScore) {
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

        ofs << seqid << "\tminiprot_scorer\tintron\t";
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
        if (introns.size() == 0 ||
            !introns[0].leftExon->initial ||
            introns[0].rightExon->score < minExonScore ||
            introns[0].score < minInitialIntronScore) {
            return;
        }
    }

    ofs << seqid << "\tminiprot_scorer\tstart_codon\t";
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
    if (introns.size() != 0 &&
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
                !introns[0].leftExon->initial ||
                introns[0].rightExon->score < minExonScore ||
                introns[0].score < minInitialIntronScore) {
                continue;
            }
        }
        ofs << seqid << "\tminiprot_scorer\tCDS\t";
        if (forward) {
            if (exons[i]->gapStart) {
                ofs << pairs[exons[i]->start].realPosition + 1 << "\t";
            } else {
                ofs << pairs[exons[i]->start].realPosition << "\t";
            }
            ofs << pairs[exons[i]->end].realPosition << "\t";
        } else {
            ofs << pairs[exons[i]->end].realPosition << "\t";
            if (exons[i]->gapStart) {
                ofs << pairs[exons[i]->start].realPosition - 1 << "\t";
            } else {
                ofs << pairs[exons[i]->start].realPosition << "\t";
            }
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
        ofs << seqid << "\tminiprot_scorer\tstop_codon\t";
        if (forward) {
            ofs << pairs[stop->position].realPosition << "\t";
            ofs << pairs[stop->position + 2].realPosition  << "\t";
        } else {
            ofs << pairs[stop->position + 2].realPosition  << "\t";
            ofs << pairs[stop->position].realPosition << "\t";
        }
        ofs << ".\t" << strand << "\t0\tprot=" << protein << ";";
        ofs << " al_score=" << stop->score << ";";
        ofs << " eScore=" << stop->exon->score << ";";

        bool proteinEnd = false;
        if (alignedProteinLength + proteinStart - 1 == proteinLength) {
            proteinEnd = true;
        }
        ofs << " proteinEnd=" << proteinEnd << ";\n";
    }
}

void Alignment::print(ostream& os) {
    for (unsigned int i = 0; i < blockLength; i++) {
        os << pairs[i].nucleotide;
        if (pairs[i].nucleotide == '<') {
            break;
        }
    }
    os << endl;
    for (unsigned int i = 0; i < blockLength; i++) {
        os << pairs[i].translatedCodon;
        if (pairs[i].translatedCodon == '<') {
            break;
        }
    }
    os << endl;
    for (unsigned int i = 0; i < blockLength; i++) {
        os << pairs[i].type;
        if (pairs[i].nucleotide == '<') {
            break;
        }
    }
    os << endl;
    for (unsigned int i = 0; i < blockLength; i++) {
        os << pairs[i].protein;
        if (pairs[i].protein == '<') {
            break;
        }
    }
    os << endl;
}

string Alignment::getSeqid() {
    return seqid;
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

    if (n == '<') {
        // Alignment end
        this->type = 'e';
        return;
    }

    if (!gapOrNT(n)) {
        cerr << "Warning: Unexpected nucleotide \"" << n << "\". " << endl;
    }

    if (n == '-') {
        this->type = 'e';
    } else if (islower(n)) {
        this->type = 'i';
    } else {
        this->type = 'e';
    }
}

double Alignment::AlignedPair::score(const ScoreMatrix * scoreMatrix,
                                     bool intronFrameshift) {
    return scoreMatrix->getScore(translatedCodon, protein, intronFrameshift);
}

Alignment::Intron::Intron() {
    scoreSet = false;
    printedLength = 0;
    realLength = -1;
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
