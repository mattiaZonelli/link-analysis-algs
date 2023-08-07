//
// Created by Mattia Zonelli on 13/07/21.
//

#ifndef IRWS_ALGORITHMS_H
#define IRWS_ALGORITHMS_H

#include <vector>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <set>
#include "CRS.h"

using std::vector;

// compute the euclidean distance between two vectors
template<typename T>
double vectors_distance(const std::vector<T> &a, const std::vector<T> &b) {
    std::vector<double> auxiliary;

    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(auxiliary),//
                   [](T element1, T element2) { return pow((element1 - element2), 2); });
    auxiliary.shrink_to_fit();

    return std::sqrt(std::accumulate(auxiliary.begin(), auxiliary.end(), 0.0));
}

// compute the hits algorithm on the given CRSmatrix
vector<float> hits(const CRSMatrix &crs_L) {
    // init authority and hub vector to column vector of 1s
    vector<float> a(crs_L.n, 1);
    vector<float> h(crs_L.n, 1);
    int i;
    for (i = 0; i < 50; i++) {
        // compute the new authority scores by doing L^T * h
        auto a_new = T_product(crs_L, h);
        // compute the new hub scores by doing L * a
        auto h_new = product(crs_L, a);

        //normalize a_new
        float sum_a = std::accumulate(a_new.begin(), a_new.end(), 0.0);
        for (int j = 0; j < a_new.size(); ++j) {
            a_new[j] = a_new[j] / sum_a;
        }
        //normalize h_new
        float sum_h = std::accumulate(h_new.begin(), h_new.end(), 0.0);
        for (int j = 0; j < h_new.size(); ++j) {
            h_new[j] = h_new[j] / sum_h;
        }

        // euclidean distance between a-a_new and h-h_new
        double error1 = vectors_distance(a, a_new);
        double error2 = vectors_distance(h, h_new);

        // update authority scores
        a = a_new;
        // update hub scores
        h = h_new;

        // if euclidean distance of a or h drop below the threashold then break
        if (error1 < (1e-10) && error2 < (1e-10)) break;

    }
    //std::cout << "HITS early exit at iter:" << i << std::endl;
    return a;

}

// find the ids of all the danglings nodes in the sparse matrix
vector<int> findDanglings(vector<vector<int>> matrix) {
    vector<int> dang;
    for (int i = 0; i < matrix.size(); ++i) {
        if (matrix[i].empty()) {
            dang.push_back(i);
        }
    }
    return dang;
}

// compute the summation for each i in Danglings do p^k[i] / n
float summation(vector<float> prob, vector<int> dang, float n) {
    float sum = 0;
    for (int i = 0; i < dang.size(); ++i) {
        sum += (prob[dang[i]] / n);
    }

    return sum;
}

// compute the page rank algorithm for the given crs matrix
vector<float> pageRank(const CRSMatrix &crs_L, const vector<vector<int>> &matrix, double d) {
    // find the ids of the danglings
    vector<int> dang = findDanglings(matrix);

    float n = matrix.size();
    // init the probability vector
    vector<float> p(n, 1 / n);

    //double d = 0.85;
    // compute the first part of the eqaution to find p^{k+1}
    double support = (1 - d) / n;
    double error = 10e-10;

    int iter;
    for (iter = 0; iter < 101 && error >= 10e-10; ++iter) {

        // product between A^t and p^k
        auto T = T_productPR(crs_L, p, matrix);
        vector<float> p_new(T.size(), 0);

        // get the summation of the probabilities of the danglings
        float sum = summation(p, dang, n);

        // complete the sum to get p^{k+1]
        for (int j = 0; j < T.size(); ++j) {
            p_new[j] = ((T[j] + sum) * d) + support;
        }

        //euclidean distance
        error = vectors_distance(p, p_new);

        //update prob distrib
        p = p_new;
/*
        if (error < (10^(-10))) {
            std::cout << "early exit at iter:" << iter << std::endl;
            break;
        }
*/
    }
    return p;
}

// for all the nodes of crs matrix return the number of inlinks
vector<float> inDegree(const CRSMatrix &A) {
    vector<float> vResults(A.m, 0);
    for (int col : A.colIndex) {
        vResults[col] += 1;
    }
    return vResults;
}

using namespace std;
// return the topK pages in decreasing order of score
set<int> topK(vector<float> &vi, int K) {
    // create the priority queue
    std::priority_queue<std::pair<float, int>> pq;
    for (int i = 0; i < vi.size(); ++i) {
        pq.push(std::make_pair(vi[i], i));
    }

    // get only the topK
    set<int> vTopK;
    for (int i = 0; i < K; ++i) {
        std::pair<float, int> top = pq.top();
        // to print the first top K pages with according score in decreasing order, uncomment following line
        std::cout << top.second << " " << top.first <<  std::endl;
        vTopK.insert(top.second);
        pq.pop();
    }
    // to print the first top K pages, uncomment the following lines
    /*for (int x : vTopK) {
        std::cout << x << " ";
    }*/
    //std::cout <<  std::endl;

    return vTopK;
}


// Function to return the intersection set of s1 and s2
set<int> intersection(set<int> s1, set<int> s2) {
    set<int> intersect;

    // Find the intersection of the two sets
    set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(intersect, intersect.begin()));

    return intersect;
}


// Function to return the Jaccard index of two sets
double jaccard_index(set<int> s1, set<int> s2) {
    // Sizes of both the sets
    double size_s1 = s1.size();
    double size_s2 = s2.size();

    // Get the intersection set
    set<int> intersect = intersection(s1, s2);

    // Size of the intersection set
    double size_in = intersect.size();

    // Calculate the Jaccard index
    // using the formula
    double jaccard_in = size_in/ (size_s1 + size_s2 - size_in);

    // Return the Jaccard index
    return jaccard_in;
}


#endif //IRWS_ALGORITHMS_H
