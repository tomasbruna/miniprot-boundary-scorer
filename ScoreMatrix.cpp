#include "ScoreMatrix.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <float.h>

using namespace std;

double ScoreMatrix::getScore(char a, char b, bool intronFrameshift) const {

    if (intronFrameshift) {
        // Penalize no matter the input if in the intron frameshift state
        return -intronFsPenalty;
    }

    a = tolower(a);
    b = tolower(b);


    if (a == '$' || b == '$' || a == '!' || b == '!' || a == '*' || b == '*') {
        return -exonFsStopPenalty;
    }

    if (a == '-' || b == '-') {
        return -gapPenalty;
    }

    if (matrix.find(a) != matrix.end()) {
        auto score = matrix.at(a).find(b);
        if (score != matrix.at(a).end()) {
            return score->second;
        }
    }

    // This should not be occuring anymore as frameshifts are handled
    // explicitly
    cerr << "Warning: score for (" << a << "," << b << ") is not defined in the "
            "matrix. Returning " << UNKNOWN_SCORE << " instead." << endl;

    return UNKNOWN_SCORE;
}

double ScoreMatrix::getMaxScore() const{
    return maxScore;
}

void ScoreMatrix::computeMaxScore() {
    maxScore = -1 * DBL_MAX;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            double score = matrix.at(columnHeaders[i]).at(columnHeaders[j]);
            if (score > maxScore) {
                maxScore = score;
            }
        }
    }
}

void ScoreMatrix::setPenalties(double gapPenalty, double intronFsPenalty,
                               double exonFsStopPenalty) {
    this->gapPenalty = gapPenalty;
    this->intronFsPenalty = intronFsPenalty;
    this->exonFsStopPenalty = exonFsStopPenalty;
}


bool ScoreMatrix::loadFromFile(string filename) {
    size = 0;
    inputStream.open(filename.c_str());
    if (!inputStream) {
        cerr << "error: Failed to open matrix file \"" << filename << "\"" << endl;
        return false;
    }

    if (!readColumnHeaders()) {
        return false;
    }

    for (int i = 0; i < size; i++) {
        if (!readRow()) {
            cerr << "error: Could not read matrix file" << endl;
            return false;
        }
    }

    inputStream.close();
    computeMaxScore();
    return true;
}

bool ScoreMatrix::readRow() {
    string line;
    if (!getline(inputStream, line)) {
        return false;
    }

    processLine(line);
    stringstream ss(line);
    char rowHeader;
    ss >> rowHeader;
    rowHeader = tolower(rowHeader);

    for (int i = 0; i < size; i++) {
        double score;
        if (!(ss >> score)) {
            return false;
        }
        matrix[rowHeader][columnHeaders[i]] = score;
    }

    return true;
}

bool ScoreMatrix::readColumnHeaders() {
    string line;
    if (!getline(inputStream, line)) {
        cerr << "error: Could not read matrix file" << endl;
        return false;
    }

    // Skip initial comments
    while (line[0] == '#') {
        if (!getline(inputStream, line)) {
            cerr << "error: Could not read matrix file" << endl;
            return false;
        }
    }

    processLine(line);
    stringstream ss(line);
    char columnHeader;
    while (ss >> columnHeader) {
        columnHeader = tolower(columnHeader);
        size++;
        columnHeaders.push_back(columnHeader);
    }
    return true;
}

void ScoreMatrix::processLine(string & line) {
    std::replace(line.begin(), line.end(), '"', ' ');
    std::replace(line.begin(), line.end(), ',', ' ');
    std::replace(line.begin(), line.end(), ';', ' ');
    std::replace(line.begin(), line.end(), '|', ' ');
}

void ScoreMatrix::print() const {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            cout << matrix.at(columnHeaders[i]).at(columnHeaders[j]) << " ";
        }
        cout << endl;
    }
}
