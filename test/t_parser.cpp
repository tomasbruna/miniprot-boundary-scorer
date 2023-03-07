#include "common.h"
#include "catch.hpp"
#include "../Parser.h"
#include <stdio.h>
#include <string>
#include <iostream>

int returnDiff(string expected, string result) {
    expected = ROOT_PATH + "/test_files/" + expected;
    return system(("diff " + expected + " " + result + " >/dev/null").c_str());
}

TEST_CASE("Test whole program with different settings") {
    Parser fileParser;
    fileParser.setWindowLegth(10);
    system(("gunzip -k " + ROOT_PATH + "/test_files/test_1.aln.gz").c_str());

    string inputFile = ROOT_PATH + "/test_files/test_1.aln";
    string output = ROOT_PATH + "/test_files/test_result";
    ScoreMatrix * scoreMatrix = new ScoreMatrix();
    scoreMatrix->loadFromFile(ROOT_PATH + "/test_files/blosum62_1.csv");
    scoreMatrix->setPenalties(4, 4, 20);
    Kernel * triangularkernel = new TriangularKernel();
    fileParser.setScoringMatrix(scoreMatrix);
    fileParser.setKernel(triangularkernel);

    SECTION("All") {
        freopen(inputFile.c_str(), "r", stdin);
        fileParser.setMinExonScore(-999999);
        fileParser.setMinInitialExonScore(-999999);
        fileParser.setMinInitialIntronScore(-999999);
        fileParser.parse(output);
        int result = returnDiff("test_parser_all.gff", output);
        CHECK(result == 0);
    }

    std::cin.clear();

    SECTION("Standard filters") {
       freopen(inputFile.c_str(), "r", stdin);
       fileParser.setMinExonScore(25);
       fileParser.setMinInitialIntronScore(0);
       fileParser.setMinInitialExonScore(25);
       fileParser.parse(output);
       int result = returnDiff("test_parser.gff", output);
       CHECK(result == 0);
    }

    delete scoreMatrix;
    delete triangularkernel;
    remove(output.c_str());
    remove((ROOT_PATH + "/test_files/test_1.aln").c_str());
}
