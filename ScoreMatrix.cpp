#include "ScoreMatrix.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <float.h>

using namespace std;

double ScoreMatrix::getScore(char a, char b) const {

    a = tolower(a);
    b = tolower(b);
    // Treat spaces and dashes as gaps
    if (a == ' ' || a == '-') {
        a = '*';
    }
    if (b == '-' || b == ' ') {
        b = '*';
    }

    if (matrix.find(a) != matrix.end()) {
        auto score = matrix.at(a).find(b);
        if (score != matrix.at(a).end()) {
            return score->second;
        }
    }

    // This is caused by frameshifts and occasional shift of how gaps are printed
    // in alignment (Spaln specific). In any case, penalty -4 is reasonable to penalize
    // both of these cases -- if it really is a frameshift, the rest of the alignment
    // also receives -4 penalty. If it is the gap case, the rest of the alignment
    // after gap is printed correctly.
    // cerr << "Warning: score for (" << a << "," << b << ") is not defined in the "
    //        "matrix. Returning " << UNKNOWN_SCORE << " instead." << endl;

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
