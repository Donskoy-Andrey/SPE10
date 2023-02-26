#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include "solver_system/bicgstab.h"
#include "solver_system/amg_ruge_stuben.h"
#include "solver_system/gauss_seidel.h"
#include "solver_system/ilduc.h"
#include "solver_system/cpr.h"
#include "solver_system/two_stage.h"

/*
---------------------------------------------------
 COO - format: ia - rows, ja - columns, a- values
 EXAMPLE:
 matrix:
    3    5    0    0
    0    6    0    0
    0    0    0    1
    0    0    0    4

 COO-format:
 ia = [0, 0, 1, 2, 3]
 ja = [0, 1, 1, 3, 3]
  a = [3, 5, 6, 1, 4]
------------------------------------------------------- 
*/


class COO {
public:

    std::vector<int> ia = {};
    std::vector<int> ja = {};
    std::vector<double> a = {};


    COO() {}

    /* Return number of rows of matrix [if matrix NxN returns N] */
    int len_mat() {
        if (ia[ia.size() - 1] > ja[ja.size() - 1]) { return ia[ia.size() - 1] + 1; }  // !!!!!!!!maybe  wrong!!!!!!!!!!
        else return ja[ja.size() - 1] + 1;
    }

    /* Insert value into the matrix */
    void insert_val(int row, int col, double value) {
        if (value == 0) { return; }
        if ((ia.empty()) ||
            ((row >= ia[ia.size() - 1]) && (col > ja[ja.size() - 1]))) {
            ia.push_back(row);
            ja.push_back(col);
            a.push_back(value);
        } else {

            const auto p = std::equal_range(ia.begin(), ia.end(), row);
            const auto index_of_first = p.first - ia.begin();
            const auto index_of_last = p.second - ia.begin();


            const auto first = next(ja.begin(), index_of_first);
            const auto last = next(ja.begin(), index_of_last);


            auto col_pos_it = upper_bound(first, last, col);
            auto pos = col_pos_it - ja.begin();
            ja.insert(col_pos_it, col);

            auto row_pos_it = next(ia.begin(), pos);
            ia.insert(row_pos_it, row);

            auto val_pos_it = next(a.begin(), pos);
            a.insert(val_pos_it, value);
        }
    }


    /* Return value by index in COO-matrix [only active elements - not zeros] */
    double operator()(int index) {
        return a[index];
    }


    double operator()(int row, int col) {
        if (row > len_mat() - 1 || row > len_mat() - 1) {
            std::cerr << "out of range";
            throw;
        } else {
            int i = 0;
            while (ia[i] != row && i <= ia.size() - 1) {
                i++;
            }
            if (i == ia.size()) { return 0; }
            while (ja[i] != col && ia[i] == row) {
                i++;
            }
            if (ia[i] != row) { return 0; }
            return a[i];
        }
    }

    /* Print matrix in coo format [ia - rows, ja - columns, a - values] */
    void print_coo() {
        std::cout << "ia = [";
        for (int i = 0; i < ia.size() - 1; ++i) {
            std::cout << ia[i] << ", ";
        }
        std::cout << ia[ia.size() - 1] << "]" << std::endl;

        std::cout << "ja = [";
        for (int i = 0; i < ja.size() - 1; ++i) {
            std::cout << ja[i] << ", ";
        }
        std::cout << ja[ja.size() - 1] << "]" << std::endl;

        std::cout << "a = [";
        for (int i = 0; i < a.size() - 1; ++i) {
            std::cout << a[i] << ", ";
        }
        std::cout << a[a.size() - 1] << "]" << std::endl;
    }


    std::vector<int> get_ia() {
        std::vector<int> ia_for_ret = ia;
        return ia_for_ret;
    }

    std::vector<int> get_ja() {
        std::vector<int> ja_for_ret = ja;
        return ja_for_ret;
    }

    std::vector<double> get_a() {
        std::vector<double> a_for_ret = a;
        return a_for_ret;
    }

void clear() {
    ia = {};
    ja = {};
    a = {};
    return;
}

/* Print matrix */
void print_matrix() {
    for (int i = 0; i < len_mat(); ++i) {
        for (int j = 0; j < len_mat(); ++j) {
            std::cout << std::setw(5) << this->operator()(i, j);
        }
        std::cout << std::endl;
    }
}

// void translate_csr(std::vector<idx_t> &iia, std::vector<idx_t> &jja, std::vector<double> &aa){  
//         int *buf = new int[ja[ja.size() - 1] + 2];
//         for (int i = 0; i < ja[ja.size() - 1] + 2; i++){
//             buf[i] = 0;
//         }
//         for (int i = 0; i < a.size(); i++)
//         {
//             jja.push_back(ja[i]);
//             aa.push_back(a[i]);

//             buf[ia[i] + 1]++;
//         }
//         std::cout << std::endl;
//         for (int i = 0; i < ja[ja.size() - 1] + 2; i++)
//         {
//             buf[i + 1] += buf[i];
//         }
//         for (int i = 0; i < ja[ja.size() - 1] + 2; i++)
//         {
//             iia.push_back(buf[i]);
//         }

//         std::ofstream file;
//         file.open("../data/A_csr.mtx");
//         for (int i = 0; i < aa.size()+1; ++i) {
//             file << iia[i] << " " << jja[i] << " " << aa[i] << std::endl;
//         }
// }

void translate_csr(std::vector<idx_t> &iia, std::vector<idx_t> &jja, std::vector<double> &aa)
{
    const int nrows = 60*220*2;
    const int ncols = a.size();

    // Step 1: Initialize the CSR vectors
    iia.resize(nrows + 1);
    jja.resize(a.size());
    aa.resize(a.size());

    // Step 2: Count the number of non-zero elements in each row
    std::vector<idx_t> nnz_row(nrows, 0);
    for (idx_t k = 0; k < a.size(); k++) {
        nnz_row[ia[k]]++;
    }

    // Step 3: Compute the prefix sum of nnz_row to obtain iia
    idx_t cumsum = 0;
    for (idx_t i = 0; i < nrows; i++) {
        iia[i] = cumsum;
        cumsum += nnz_row[i];
    }
    iia[nrows] = a.size();

    // Step 4: Fill in jja and aa
    std::vector<idx_t> next_row_idx(nrows, 0);
    for (idx_t k = 0; k < a.size(); k++) {
        idx_t row = ia[k];
        idx_t idx = iia[row] + next_row_idx[row];
        jja[idx] = ja[k];
        aa[idx] = a[k];
        next_row_idx[row]++;
    }
    #if SAVE_CSR_MATRIX
        std::ofstream file;
        file.open("../data/A_csr.mtx");
        for (int i = 0; i < iia.size()+1; ++i) {
            file << iia[i] << " ";
        }
        file << std::endl;
        for (int i = 0; i < jja.size()+1; ++i) {
            file << jja[i] << " ";
        }
        file << std::endl;
        for (int i = 0; i < aa.size()+1; ++i) {
            file << aa[i] << " ";
        }
    #endif
}


/* Save matrix in COO format to ["../data/A.mtx" = default]*/
void write_to_file(std::string filename = "../data/A.mtx") {
    std::ofstream file;
    file.open(filename);
    for (int i = 0; i < ia.size()+1; ++i) {
        file << ia[i] << " " << ja[i] << " " << a[i] << std::endl;
    }
}

};
