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

int Alignment::parse(istream& inputStream, string headerLine, bool forward) {
    clear();
    this->forward = forward;
    vector<string> blockLines(BLOCK_ITEMS_CNT);
    string line;

    // Read header
    int status = parseHeader(headerLine);
    if (status != READ_SUCCESS) {
        cerr << "error: Invalid alignment header " << endl;
        return status;
    }

    // Get to the alignment itself
    getline(inputStream, line);
    while (line.substr(0, 9) != "ALIGNMENT") {
        if (!getline(inputStream, line)) {
            cerr << "error: Alignment is missing after header " << endl;
            return FORMAT_FAIL;
        }
    }
    if (!getline(inputStream, line) || line != "") {
        cerr << "error: Empty line expected after ALIGNMENT keyword" << endl;
        return FORMAT_FAIL;
    }

    // Load the alignment into array
    for (i = 0; i < BLOCK_ITEMS_CNT; i++) {
        if (getline(inputStream, line) && !line.empty() && line.size() > BLOCK_OFFSET) {
            blockLines[i] = line;
            // Get alignment start number
            if (i == 1) {
                dnaStart = atoi(line.substr(0, BLOCK_OFFSET).c_str());
                if (dnaStart <= 0) {
                    cerr << "error: Could not read dna alignment start position" << endl;
                    return FORMAT_FAIL;
                }
                realPositionCounter = dnaStart;
            } else if (i == 2) {
                proteinStart = atoi(line.substr(0, BLOCK_OFFSET).c_str());
                if (proteinStart <= 0) {
                    cerr << "error: Could not read dna alignment start position" << endl;
                    return FORMAT_FAIL;
                }
            }
        } else {
            printLineError();
            return FORMAT_FAIL;
        }
    }

    // Trim the alignment strings for easier parsing.
    // Any symbol in DNA other than space before pipe signifies the true end of
    // alignment. Looking for the first whitespace in DNA does not work as there
    // might be gaps inside introns
    int pipePosition = blockLines[1].find("|");
    blockLength = blockLines[1].find_last_not_of(" ", pipePosition - 1) - BLOCK_OFFSET  + 1;
    for (i = 0; i < BLOCK_ITEMS_CNT; i++) {
        blockLines[i] = blockLines[i].substr(BLOCK_OFFSET, blockLength);
        if (blockLines[i].size() != blockLength) {
            printLineError();
            return FORMAT_FAIL;
        }
    }

    // Unify gaps
    replace(blockLines[1].begin(), blockLines[1].end(), ' ', '-');
    // Start actual parsing
    parseBlock(blockLines);

    // Potential single-exon gene
    if (exons.size() == 1) {
        exons.back()->initial = false;
    }

    return READ_SUCCESS;
}

void Alignment::printLineError() {
    cerr << "warning: error in alignment " << gene << "-" << protein;
    cerr << ": corrupted alignment - wrong line length.";
    cerr << " The rest of this alignment is skipped." << endl;
}

int Alignment::parseHeader(string headerLine) {
    bool geneParsed = false;
    for (unsigned int i = 0; i < headerLine.size(); i++) {
        if (headerLine[i] == '>' || headerLine[i] == '<') {
            if (!geneParsed) {
                int length = headerLine.find(" ", i + 1) - i - 1;
                gene = headerLine.substr(i + 1, length);
                geneParsed = true;
            } else {
                int length = headerLine.find(" ", i + 1) - i - 1;
                protein = headerLine.substr(i + 1, length);
                return READ_SUCCESS;
            }
        }
    }
    return FORMAT_FAIL;
}

void Alignment::parseBlock(const vector<string>& lines) {
    // Parse individual pairs
    for (unsigned i = 0; i < lines[0].size(); i++) {
        AlignedPair pair(lines[0][i], lines[1][i], lines[2][i], insideIntron);

        checkForIntron(pair);
        checkForStart(pair);
        checkForStop(pair);

        if (pair.nucleotide != '-' && pair.nucleotide != ' ') {
            if (forward) {
                pair.realPosition = realPositionCounter++;
            } else {
                pair.realPosition = realPositionCounter--;
            }
        } else {
            // Don't increment the counter
            if (forward) {
                pair.realPosition = realPositionCounter - 1;
            } else {
                pair.realPosition = realPositionCounter + 1;
            }
        }

        // Reuse space if possible
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

void Alignment::assignCodonPhases() {
    for (unsigned int i = 0; i < blockLength; i++) {
        if (pairs[i].translatedCodon == ' ' && pairs[i].type == 'e') {
            if (i == 0) {
                pairs[i].translatedCodon = '1';
            } else if (i == blockLength - 1) {
                pairs[i].translatedCodon = '3';
            } else if (gapOrAA(pairs[i + 1].translatedCodon)) {
                pairs[i].translatedCodon = '1';
            } else if (gapOrAA(pairs[i - 1].translatedCodon)) {
                pairs[i].translatedCodon = '3';
            } else if (pairs[i + 1].type == 'i') {
                pairs[i].translatedCodon = '1';
            } else if (pairs[i - 1].type == 'i') {
                pairs[i].translatedCodon = '3';
            }
        }

        if (pairs[i].protein == ' ' && pairs[i].type == 'e') {
            if (i == 0) {
                pairs[i].protein = '1';
            } else if (i == blockLength - 1) {
                pairs[i].protein = '3';
            } else if (gapOrAA(pairs[i + 1].protein)) {
                pairs[i].protein = '1';
            } else if (gapOrAA(pairs[i - 1].protein)) {
                pairs[i].protein = '3';
            } else if (pairs[i + 1].type == 'i') {
                pairs[i].protein = '1';
            } else if (pairs[i - 1].type == 'i') {
                pairs[i].protein = '3';
            }
        }
    }
}

void Alignment::checkForIntron(AlignedPair& pair) {
    // Some alignments start with a gap or intron, do not create initial exon
    // in such cases
    int alignmentPosition;
    if (forward) {
        alignmentPosition = realPositionCounter - dnaStart + 1;
    } else {
        alignmentPosition = dnaStart - realPositionCounter + 1;
    }
    if (alignmentPosition == 1 && pair.type == 'e' && pair.nucleotide != '-') {
        exons.push_back(new Exon(index));
    }

    // Alignment end
    if (index == (int) blockLength - 1 && pair.type == 'e') {
        exons.back()->end = index;
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
        exons.push_back(new Exon(index));

        // Make the decision about exon phase based on
        // how the preceeding exon was split.
        if (introns.back().start != 0) {
            if (gapOrAA(pairs[introns.back().start - 1].protein)) {
                exons.back()->phase = 1;
            } else if (introns.back().start - 2 < 0) {
                // Reached start of alignment
                exons.back()->phase = 2;
            } else {
                int position = introns.back().start - 2;
                if (pairs[position].type != 'e') {
                    // Check if still in exon (some exons are just 1 nt long).
                    // If not, take last nt from next exon upstream
                    position = introns[introns.size() - 2].start - 1;
                }
                if (position >= 0 && gapOrAA(pairs[position].protein)) {
                    exons.back()->phase = 0;
                } else {
                    exons.back()->phase = 2;
                }
            }
        }

        // Frameshift, remove false intron
        if (index - introns.back().start < 3) {
            introns.pop_back();
        } else {
            introns.back().end = index - 1;
            introns.back().acceptor[0] = pairs[index - 2].nucleotide;
            introns.back().acceptor[1] = pairs[index - 1].nucleotide;
            if (introns.back().start != 0 && introns.back().gap == false) {
                introns.back().complete = true;
            }
            introns.back().rightExon = exons.back();
        }
    } else if (insideIntron && pair.nucleotide == '-') {
        // Gap (AA aligned) inside introns, do not report these introns
        introns.back().gap = true;
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

        ofs << gene << "\tSpaln_scorer\tIntron\t";
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

    ofs << gene << "\tSpaln_scorer\tstart_codon\t";
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
        ofs << gene << "\tSpaln_scorer\tCDS\t";
        if (forward) {
            ofs << pairs[exons[i]->start].realPosition << "\t";
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
        ofs << gene << "\tSpaln_scorer\tstop_codon\t";
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
        os << pairs[i].translatedCodon;
    }
    os << endl;
    for (unsigned int i = 0; i < blockLength; i++) {
        os << pairs[i].nucleotide;
    }
    os << endl;
    for (unsigned int i = 0; i < blockLength; i++) {
        os << pairs[i].protein;
    }
    os << endl;
    for (unsigned int i = 0; i < blockLength; i++) {
        os << pairs[i].type;
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

Alignment::AlignedPair::AlignedPair(char tc, char n, char p, bool insideIntron) :
nucleotide(n),
translatedCodon(tc),
protein(p) {
    if (islower(n) || (insideIntron && n == '-')) {
        this->type = 'i';
    } else {
        this->type = 'e';
    }
    // Assign an amino acid to a stop codon, so it is treated as a gap
    if (translatedCodon == '*') {
        translatedCodon = 'A';
    }

    // Fix strange J amino acid which is in fact S
    if (translatedCodon == 'J') {
        translatedCodon = 'S';
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

Alignment::Exon::Exon(int start) {
    this->start = start;
    phase = 0;
    scoreSet = false;
    initial = false;
}

Alignment::Codon::Codon(int position, Exon * exon) {
    this->position = position;
    this->exon = exon;
}
