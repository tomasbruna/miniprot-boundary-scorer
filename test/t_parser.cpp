

#include "common.h"
#include "catch.hpp"
//#include "catch_amalgamated.hpp"
#include "../Parser.h"
#include <stdio.h>
#include <string>
#include <iostream>

// system(diff) in this testing case is system dependent
// this is ok for testing purposes

int returnDiff(string expected, string result) {
    expected = ROOT_PATH + "/test_files/" + expected;
    return system(("diff " + expected + " " + result + " >/dev/null").c_str());
}

TEST_CASE("Test whole program with different settings") {
    Parser fileParser;
    fileParser.setWindowLegth(10);
    system(("gunzip -k " + ROOT_PATH + "/test_files/test_1.ali.gz").c_str());

    string inputFile = ROOT_PATH + "/test_files/test_1.ali";
    string output = ROOT_PATH + "/test_files/test_result";
    ScoreMatrix * scoreMatrix = new ScoreMatrix();
    scoreMatrix->loadFromFile(ROOT_PATH + "/test_files/blosum62_1.csv");
    Kernel * triangularkernel = new TriangularKernel();
    fileParser.setScoringMatrix(scoreMatrix);
    fileParser.setKernel(triangularkernel);

    SECTION("All") {
        freopen(inputFile.c_str(), "r", stdin);
        fileParser.setMinExonScore(-999999);
        fileParser.setMinInitialExonScore(-999999);
        fileParser.setMinInitialIntronScore(-999999);
        fileParser.parse(output);
        int result = returnDiff("test_parser.gff", output);
        CHECK(result == 0);
    }

    std::cin.clear();

   // SECTION("Standard filters") {
     //   freopen(inputFile.c_str(), "r", stdin);
      //  fileParser.setMinExonScore(25);
      //  fileParser.setMinInitialIntronScore(0.1);
     //   fileParser.setMinInitialExonScore(0);
     //   fileParser.parse(output);
     //   int result = returnDiff("test_parser_eScore_25.gff", output);
     //   CHECK(result == 0);
   // }

    delete scoreMatrix;
    delete triangularkernel;
    //remove(output.c_str());
    remove((ROOT_PATH + "/test_files/test_1.ali").c_str());
}
