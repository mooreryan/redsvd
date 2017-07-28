/*
 *  Copyright (c) 2010 Daisuke Okanohara
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 */

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "redsvdFile.hpp"
#include "redsvd.hpp"
#include "redsvdIncr.hpp"

using namespace std;
using namespace Eigen;

/* magnitude of vector */
double mag(VectorXf v)
{
  return sqrt(v.array().square().sum());
}

double cosine_similarity(VectorXf v1, VectorXf v2)
{
  return v1.dot(v2) / (mag(v1) * mag(v2));
}

/* TODO other options may be better eg angular distance */
double cosine_dissimilarity(VectorXf v1, VectorXf v2)
{
  return 1 - abs(cosine_similarity(v1, v2));
}

namespace REDSVD{

  namespace {

    void writeAVADis_(const string& fn, const MatrixXf& M){
      cout << "write " << fn << endl;
      FILE* outfp = fopen(fn.c_str(), "wb");
      if (outfp == NULL){
        throw string("cannot open ") + fn;
      }

      int num_rows = M.rows();

      for (int r1 = 0; r1 < num_rows - 1; ++r1) {
        if (r1 % 100 == 0) {
          fprintf(stderr,
                  "Writing row %d of %d\r",
                  r1,
                  num_rows);
        }

        for (int r2 = r1 + 1; r2 < num_rows; ++r2) {
          fprintf(outfp, "%d %d %+f\n", r1, r2,
                  cosine_dissimilarity(M.row(r1), M.row(r2)));
        }
      }

      fclose(outfp);

    }

    /* Rows will be rows from M1 and columns will be the dissimilarity
       scores of this row with each row in M2. So this 2nd row, 3rd
       column would be the dissimilarity between M1.row(1) and
       M2.row(2). */
    void writeABDis_(const string& fn,
                     const MatrixXf& M1,
                     const MatrixXf& M2){

      cout << "write " << fn << endl;
      FILE* outfp = fopen(fn.c_str(), "wb");
      if (outfp == NULL){
        throw string("cannot open ") + fn;
      }

      double cos_dis = 0;
      vector<double> dists;
      vector<double>::iterator iter;

      for (int i1 = 0; i1 < M1.rows(); ++i1) {
        for (int i2 = 0; i2 < M2.rows(); ++i2) {
          cos_dis = cosine_dissimilarity(M1.row(i1), M2.row(i2));

          dists.push_back(cos_dis);
        }

        iter = dists.begin();
        fprintf(outfp, "%+f", *iter);
        for (iter = dists.begin() + 1; iter < dists.end(); ++iter) {
          fprintf(outfp, " %+f", *iter);
        }
        fprintf(outfp, "\n");
        dists.clear();
      }

      fclose(outfp);
    }


    void writeMatrix_(const string& fn, const MatrixXf& M){
      cout << "write " << fn << endl;
      FILE* outfp = fopen(fn.c_str(), "wb");
      if (outfp == NULL){
        throw string("cannot open ") + fn;
      }

      for (int i = 0; i < M.rows(); ++i){
        /* print first element of the line */
        fprintf(outfp, "%+f",  M(i, 0));

        for (int j = 1; j < M.cols(); ++j){
          fprintf(outfp, " %+f",  M(i, j));
        }
        fprintf(outfp, "\n");
      }

      fclose(outfp);
    }

    void writeVector_(const string& fn, const VectorXf& V){
      cout << "write " << fn << endl;
      FILE* outfp = fopen(fn.c_str(), "wb");
      if (outfp == NULL){
        throw string("cannot open ") + fn;
      }

      for (int i = 0; i < V.rows(); ++i){
        fprintf(outfp, "%+f\n", V(i));
      }

      fclose(outfp);
    }

    void readLine(const string& line,
                  fv_t& fv){
      istringstream is(line);

      int id;
      char sep;
      float val;
      while (is >> id >> sep >> val){
        fv.push_back(make_pair(id, val));
      }
      sort(fv.begin(), fv.end());
      fv.erase(unique(fv.begin(), fv.end()), fv.end());
    }

  }



  void readMatrix(const std::string& fn, SMatrixXf& A){
    vector<fv_t> fvs;
    ifstream ifs(fn.c_str());
    if (!ifs){
      throw string("failed to open") + fn;
    }

    for (string line; getline(ifs, line); ){
      fv_t fv;
      readLine(line, fv);
      //if (fv.size() == 0) continue;
      fvs.push_back(fv);
    }
    Util::convertFV2Mat(fvs, A);
  }

  void readMatrix(const std::string& fn, MatrixXf& A){
    ifstream ifs(fn.c_str());
    if (!ifs){
      throw string("failed to open " ) + fn;
    }

    vector< vector<float> > vs;
    for (string line; getline(ifs, line); ){
      istringstream is(line);
      vector<float> v;
      float val;
      while (is >> val){
        v.push_back(val);
      }
      vs.push_back(v);
    }

    size_t rowN = vs.size();
    if (rowN == 0) return;
    size_t colN = vs[0].size();
    A.resize(rowN, colN);

    for (size_t i = 0; i < rowN; ++i){
      if (colN != vs[i].size()){
        cerr << "warning: " << i+1 << "-th row has " << vs[i].size() << " entries. "
             << colN << " entries are expected" << endl;
      }
      size_t colNmin = min(colN, vs[i].size());
      for (size_t j = 0; j < colNmin; ++j){
        A(i, j) = vs[i][j];
      }
    }
  }

  void writeMatrix(const string& fn, const REDSVD::RedSVD& A){
    writeMatrix_(fn + ".U", A.matrixU());
    writeVector_(fn + ".S", A.singularValues());
    writeMatrix_(fn + ".V", A.matrixV());

    /* And these are the terms in latent space in the context of
       lsa.rb */
    MatrixXf US = A.matrixU() * A.singularValues().asDiagonal();

    /* In the context of lsa.rb, these are the documents in latent
       space */
    MatrixXf VS = A.matrixV() * A.singularValues().asDiagonal();

    writeMatrix_(fn + ".US", US);
    writeMatrix_(fn + ".VS", VS);

    /* This can get expensive if there are a lot of terms, which there
       will be. So only calculate it on the docs and the docs vs
       terms. */

    if (US.rows() < VS.rows()) {
      writeAVADis_(fn + ".US.dis", US);
      writeABDis_(fn + ".US_to_VS.dis", US, VS);
    } else if (VS.rows() < US.rows()) {
      writeAVADis_(fn + ".VS.dis", VS);
      writeABDis_(fn + ".VS_to_US.dis", VS, US);
    } else { /* if they are the same, just print both */
      writeAVADis_(fn + ".US.dis", US);
      writeABDis_(fn + ".US_to_VS.dis", US, VS);

      writeAVADis_(fn + ".VS.dis", VS);
      writeABDis_(fn + ".VS_to_US.dis", VS, US);
    }
  }

  void writeMatrix(const string& fn, const REDSVD::RedSVDIncr& A){
    writeMatrix_(fn + ".U", A.matrixU());
    writeVector_(fn + ".S", A.singularValues());
    writeMatrix_(fn + ".V", A.matrixV());
  }


  void writeMatrix(const string& fn, const REDSVD::RedPCA& A){
    writeMatrix_(fn + ".pc",    A.principalComponents());
    writeMatrix_(fn + ".score", A.scores());
  }

  void writeMatrix(const string& fn, const REDSVD::RedSymEigen& A){
    writeMatrix_(fn + ".evec", A.eigenVectors());
    writeVector_(fn + ".eval", A.eigenValues());
  }

}
