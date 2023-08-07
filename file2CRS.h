//
// Created by Mattia Zonelli on 07/07/21.
//

#ifndef IRWS_FILE2CRS_H
#define IRWS_FILE2CRS_H

#include "CRS.h"

#include <string>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstdlib>
#include <ctime>

using std::ifstream;
using std::cerr;
using std::string;

vector<vector<int>> read(const string &filename) {
    string line;
    //open the file
    ifstream input_file(filename);
    if (!input_file.is_open()) {
        cerr << "Could not open the file -> " << filename << "" << std::endl;
        return std::vector<std::vector<int>>(10);
    }

    // skip first three lines
    for (int i = 0; i < 3; ++i) {
        getline(input_file, line);
    }

    //std::cout << line << std::endl;
    //get the number of nodes
    vector<string> cut(5);
    std::stringstream str(line);
    str >> cut[0];
    str >> cut[1];
    str >> cut[2]; // questo è il numero di nodi

    int n_nodes = std::stoi(cut[2]);
    //int n_edges = std::stoi(cut[4]);

    // skip another line of the file
    getline(input_file, line);
    // check the correct number of nodes
    while (getline(input_file, line)) {
        string strfromId, strtoId;
        std::stringstream str_tmp(line);
        str_tmp >> strfromId;
        str_tmp >> strtoId;
        int fromId = std::stoi(strfromId);
        int toId = std::stoi(strtoId);
        if (fromId > n_nodes) n_nodes = fromId + 1;
        if (toId > n_nodes) n_nodes = toId + 1;
    }
    // go back to the first line of the file
    input_file.clear();
    input_file.seekg(0);

    // skip first 4 lines
    for (int i = 0; i < 4; ++i) {
        getline(input_file, line);
    }

    // create the adjacency matrix as a vector of vector
    // in position i there is a vector with the ids of the pointed nodes by node i
    std::vector<std::vector<int>> matrix(n_nodes);
    while (getline(input_file, line)) {
        string strfromId, strtoId;
        std::stringstream str_tmp(line);
        str_tmp >> strfromId;
        str_tmp >> strtoId;
        int fromId = std::stoi(strfromId);
        int toId = std::stoi(strtoId);
        //matrix[fromId-1].push_back(toId);
        matrix[fromId].push_back(toId);
    }
    input_file.close();
    return matrix;
}

CRSMatrix toCRS(const vector<vector<int>> &matrix) {

    int matrix_len = matrix.size();
    int nz = 0;
    for (int i = 0; i < matrix_len; ++i) {
        nz+=matrix[i].size();
    }

    // init the three main vector of crs representation
    vector<float> row_ptr(matrix_len+1, 0.0);
    vector<float> val(nz, 1.0), col_ind(nz, 0.0);

    float count = 0;
    int i;
    // scorriamo all the matrix
    for (i = 0; i < matrix_len; ++i) {
        // save the pointer to the first cell of column index in that row
        row_ptr[i] = count;
        if (!matrix[i].empty()) {
            // for all the element pointed by i
            for (int j = 0; j < matrix[i].size(); ++j) {
                count++;
                int tmp = int (count-1);
                // save the column index of the pointed node
                col_ind[tmp] = matrix[i][j];
            }
        }
    }
    // save one more  index of the actual number of rows)
    row_ptr[i] = count;

    CRSMatrix res{
            (int) matrix_len, (int) matrix_len, (int) count,
            val, col_ind, row_ptr
    };
    return res;
}

/*
CRSMatrix sparseMatrices(const vector<vector<int>> &matrix) {
    
    FILE *fstream_col, *fstream_val, *fstream_row;

    fstream_col = fopen("./mmapped_file_col", "w+");
    fstream_val = fopen("./mmapped_file_val", "w+");
    fstream_row = fopen("./mmapped_file_row", "w+");

    float count = 0, count_row = 0;
    for (const auto &fromId : matrix) {
        fwrite(&count, sizeof(float), 1, fstream_row);
        if (!fromId.empty()) {
            for (auto toId : fromId) {
                count++;
                float f_toId = toId;
                fwrite(&f_toId, sizeof(float), 1, fstream_col);

                float uno = 1.0;
                fwrite(&uno, sizeof(float), 1, fstream_val);
            }
        }

        count_row++;

    }
    fwrite(&count, sizeof(float), 1, fstream_row);
    

    fclose(fstream_col);
    fclose(fstream_val);
    fclose(fstream_row);

    float *mmap_region_col;
    int fd_col;
    fd_col = open("./mmapped_file_col", O_RDONLY);
    mmap_region_col = (float *) mmap(nullptr, count * sizeof(float), PROT_READ, MAP_SHARED, fd_col, 0);
    if (mmap_region_col == MAP_FAILED) {
        close(fd_col);
        printf("Error mmapping the file");
        exit(1);
    }
    close(fd_col);

    float *mmap_region_val;
    int fd_val;
    fd_val = open("./mmapped_file_val", O_RDONLY);
    mmap_region_val = (float *) mmap(nullptr, count * sizeof(float), PROT_READ, MAP_SHARED, fd_val, 0);
    if (mmap_region_val == MAP_FAILED) {
        close(fd_val);
        printf("Error mmapping the file");
        exit(1);
    }
    close(fd_val);

    float *mmap_region_row;
    int fd_row;
    fd_row = open("./mmapped_file_row", O_RDONLY);
    mmap_region_row = (float *) mmap(nullptr, count_row * sizeof(float), PROT_READ, MAP_SHARED, fd_row, 0);
    if (mmap_region_row == MAP_FAILED) {
        close(fd_row);
        printf("Error mmapping the file");
        exit(1);
    }
    close(fd_row);


    CRSMatrix res{
            (int) count_row, (int) count_row, (int) count,
            mmap_region_val, mmap_region_col, mmap_region_row
    };


    return res;
}
*/

#endif //IRWS_FILE2CRS_H
