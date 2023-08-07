#include "file2CRS.h"
#include "algorithms.h"
#include <fstream>


int main(int argc, char *argv[]) {

    //std::string string1("hitsTest");
    //auto matrix = read("/Users/mattia/CLionProjects/irws/data/" + string1 + ".txt");

    auto matrix = read("./data/hitsTest.txt");

    CRSMatrix crs_L = toCRS(matrix);

    auto res3 = hits(crs_L);
    //auto res2 = pageRank(crs_L, matrix, 0.85);
    //auto res1 = inDegree(crs_L);

    auto top = topK(res3, 8);
/*
    ofstream myfile;
    myfile.open(string1 + "-inDeg-PR.txt");
    myfile << "K" << "\t" << "jacc_coeff" << "\n";

    for (int i = 10; i <= 100; i += 10) {
        auto top1 = topK(res1, i);
        auto top2 = topK(res2, i);
        auto jacc_coeff = jaccard_index(top1, top2);
        //std::cout << "Jaccard Coefficient : " << jacc_coeff  << " of top: " << i << std::endl;
        myfile << i << "\t" << jacc_coeff << "\n";
    }
    myfile.close();

    ofstream myfile1;
    myfile1.open(string1 + "-PR-hits.txt");
    myfile1 << "K" << "\t" << "jacc_coeff" << "\n";
    for (int i = 10; i <= 100; i += 10) {
        auto top1 = topK(res2, i);
        auto top2 = topK(res3, i);
        auto jacc_coeff = jaccard_index(top1, top2);
        //std::cout << "Jaccard Coefficient : " << jacc_coeff  << " of top: " << i << std::endl;
        myfile1 << i << "\t" << jacc_coeff << "\n";
    }
    myfile1.close();

    ofstream myfile2;
    myfile2.open(string1+"-hits-inDeg.txt");
    myfile2 << "K" << "\t" << "jacc_coeff" << "\n";
    for (int i = 10; i <= 100; i+=10) {
        auto top1 = topK(res1, i);
        auto top2 = topK(res3, i);
        auto jacc_coeff = jaccard_index(top1, top2);
        //std::cout << "Jaccard Coefficient : " << jacc_coeff  << " of top: " << i << std::endl;
        myfile2 << i << "\t" << jacc_coeff << "\n";
    }
    myfile2.close();
*/
}