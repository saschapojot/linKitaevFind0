//
// Created by tusa on 5/1/21.
//

#ifndef LINKITAEVFIND0_SOLPRODUCER_HPP
#define LINKITAEVFIND0_SOLPRODUCER_HPP
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <map>
#include "consts.hpp"
#include <complex>
#include <fstream>
#include <cmath>
#include <utility>
#include <algorithm>
#include <memory>
#include <cmath>
#include <numeric>
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
class producer {
public:
    producer(const CONSTS &);
/*
 * coefs before the quench
 **/
    double b0(const double &k);

    double c0(const double&k);
    /*
     * coefs after the quench
     * */

    double b1(const double&k);

    double c1(const double&k);
    /*
     * init eigenvector at momentum k
     * */

    Eigen::Vector2cd phi0(const double&k);
/*
 *
 * time evolution*/
    Eigen::Matrix2cd U1(const double& t, const double &k);


/*
 * return one z value
 * */

    std::complex<double> calZ(const double &s, const double &k);

    std::vector<std::pair<int, std::complex<double>>> kAndSortZ(const double &s);

    void print123(const double&s);



public:
    CONSTS CON;
   //std::vector<std::complex<double>>zVals;

};

#endif //LINKITAEVFIND0_SOLPRODUCER_HPP
