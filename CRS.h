#include <climits>
//
// Created by Mattia Zonelli on 09/07/21.
//

#ifndef IRWS_CRS_H
#define IRWS_CRS_H

#include <vector>
#include <numeric>
#include <iostream>

using std::vector;

struct CRSMatrix{
    int n; // number of rows
    int m; // number of columns
    int nz; // number of non-zero elements
/*
    float* val; // non-zero elements
    float* colIndex; // column indices
    float* rowPtr; // row ptr
*/
    vector<float> val, // non-zero elements
    colIndex, // column indices
    rowPtr; // row pointers
};


/*
 * return the value in position [row, col]
 */
float CRS_getVal(int row, int col, const CRSMatrix& input){
    int ptr = (int)input.rowPtr[row];
    int gap = (int)input.rowPtr[row+1] - ptr; // non zero elems of desired row
    for (int i = 0; i < gap && input.colIndex[ptr+i] <= col ; ++i) {
        if (input.colIndex[ptr+i] == col){
            return input.val[ptr+i];
        }
    }
    return 0;
}

/*
 * return the value in position [row, col] considering the Trasnpose Matrix
 */
float CRS_getValTranspose(int row, int col, const CRSMatrix& input){
    return CRS_getVal(col, row, input);
}

/*
 * return the result of the product between the CRSmatrix and a column vector
 */
vector<float> product(const  CRSMatrix& A, vector<float> vi){
    // init the result vector to 0
    vector<float>vResults(vi.size(),0);
    // we need to visit all the cell of the rowPtr vector
    for (int i = 0; i < A.n; ++i) {
        // check if row i is non-zero
        if (A.rowPtr[i]!=A.rowPtr[i+1]){

            int start = A.rowPtr[i];
            int end = A.rowPtr[i+1];
            float result_acc = 0; // where we store the sum of the actual row
            // we multiply each element of the row by according element of the column vector and then we accumulate
            // the results
            for (int j = start; j < end; ++j) {
                result_acc += A.val[j] * vi[A.colIndex[j]];
            }
            vResults[i] = result_acc;
        }
    }
    return vResults;
}

/*
 * return the result of the product between the CRSmatrix consider as its transpose and a column vector
 */
vector<float> T_product(const  CRSMatrix& A, vector<float> vi){
    // init the result vector to 0
    vector<float>vResults(vi.size(),0);
    // visit matrix row by row
    // multiply each value in row i by value in pos i of column vector (so same row for both)
    // accumulate in the result vector the result of multiplication in pos of the current column of the row in matrix
    // next row ( for both matrix e result vector), repeat
    for (int i = 0; i < A.m; ++i) {
        int start = A.rowPtr[i];
        int end = A.rowPtr[i+1];
        for (int j = start; j < end; ++j) {
            vResults[A.colIndex[j]] += (A.val[j] * vi[i]);
        }
    }
    return vResults;
}

/*
 * return the result of the product between the CRSmatrix A transpose and a column vector
 * where A is the Adjacency matrix.
 */
vector<float> T_productPR(const  CRSMatrix& A, vector<float> vi, vector<vector<int>> matrix){
    vector<float>vResults(vi.size(),0);
    for (int i = 0; i < A.m; ++i) {
        if (A.rowPtr[i] != A.rowPtr[i+1]){
            int start = A.rowPtr[i];
            int end = A.rowPtr[i+1];
            for (int j = start; j < end; ++j) {
                // instead of multiply the actual values of the matrix, we use 1/(out-degree of node i)
                vResults[A.colIndex[j]] += (1.0/matrix[i].size() * vi[i]);
            }
        }
    }
    return vResults;
}

// print the CSR as a sparse matrix
void printSparseMatrix(const CRSMatrix &input) {
    for (int j = 0; j < input.m * 3 + 2; ++j) {
        std::cout << "-";
    }
    std::cout << "\n";
    for (int i = 0; i < input.n; i++) {
        std::cout << "|";
        for (int j = 0; j < input.m; j++) {
            int cell = CRS_getVal(i, j, input);
            if (cell < 10 && cell >= 0)
                std::cout << " " << cell << " ";
            else
                std::cout << cell << " ";
        }
        std::cout << "|\n";

    }
    for (int j = 0; j < input.m * 3 + 2; ++j) {
        std::cout << "-";
    }
    std::cout << std::endl;
}

// print the transpose  as a sparse matrix of the CRSmatrix
void printTranspose(const CRSMatrix &input) {
    for (int j = 0; j < input.m * 3 + 2; ++j) {
        std::cout << "-";
    }
    std::cout << "\n";
    for (int i = 0; i < input.n; i++) {
        std::cout << "|";
        for (int j = 0; j < input.m; j++) {
            int cell = CRS_getValTranspose(i, j, input);
            if (cell < 10 && cell >= 0)
                std::cout << " " << cell << " ";
            else
                std::cout << cell << " ";
        }
        std::cout << "|\n";

    }
    for (int j = 0; j < input.m * 3 + 2; ++j) {
        std::cout << "-";
    }
}

#endif //IRWS_CRS_H
