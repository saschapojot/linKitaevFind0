//
// Created by tusa on 5/1/21.
//

#include "solProducer.hpp"
producer::producer(const CONSTS &conTmp) {
    //params
    this->CON=conTmp;
}

double producer::b0(const double&k) {
    return this->CON.d0 * std::sin(k);

}
double producer::c0(const double&k) {
    return this->CON.t0 * std::cos(k) + this->CON.mu0 / 2;

}

Eigen::Vector2cd producer::phi0(const double&k) {
    double b0Val = this->b0(k);
    double c0Val = this->c0(k);
    double sum2Tmp = b0Val*b0Val + c0Val*c0Val;
    double nn = 2 * sum2Tmp - 2 * c0Val * std::sqrt(sum2Tmp);
    Eigen::Vector2cd rst;
    if (nn > this->CON.tol) {
        std::complex<double> numer1Tmp = std::complex<double>{0, b0Val};
        std::complex<double> numer2Tmp = std::complex<double>{c0Val - std::sqrt(sum2Tmp), 0};
        double denom = std::sqrt(nn);
        rst[0] = numer1Tmp / denom;
        rst[1] = numer2Tmp / denom;
        return rst;
    } else {
        if (c0Val > 0) {
            rst[0] = std::complex<double>(0, 1);
            rst[1] = std::complex<double>(0, 0);
            return rst;

        } else {
            rst[0] = std::complex<double>(0, 0);
            rst[1] = std::complex<double>(1, 0);
            return rst;


        }
    }


}

double producer::b1(const double&k) {
    return this->CON.d1 * std::sin(k);


}

double producer::c1(const double&k) {
    return this->CON.t1 * std::cos(k) + this->CON.mu1 / 2;
}

Eigen::Matrix2cd producer::U1(const double&s, const double&k) {
    Eigen::Matrix2cd I, s2, s3;
    //complex unit
    std::complex<double> i=std::complex<double>{0,1};
    //intialize id matrix I
    I(0, 0) = 1;
    I(0, 1) = 0;
    I(1, 0) = 0;
    I(1, 1) = 1;
    //initialize Pauli matrix s2
    s2(0, 0) = 0;
    s2(0, 1) = -i;
    s2(1, 0) = i;
    s2(1, 1) = 0;
    //initialize Pauli matrix s3
    s3(0, 0) = 1;
    s3(0, 1) = 0;
    s3(1, 0) = 0;
    s3(1, 1) = -1;
    double b1Val = this->b1(k);
    double c1Val = this->c1(k);

    double b12c12tmph = std::sqrt(b1Val * b1Val + c1Val * c1Val);
    Eigen::Matrix2cd V1 = std::cos(s * b12c12tmph) * I +
                          i * (b1Val / b12c12tmph * s2 + c1Val / b12c12tmph * s3) *
                          std::sin(s * b12c12tmph
                          );
    Eigen::Matrix2cd U1 = std::exp(i * this->CON.mu1 * s / 2.0) * V1;
    return U1;


}

std::complex<double> producer::calZ(const double &s, const double &k) {
    Eigen::Vector2cd phi0k = this->phi0(k);
    Eigen::Matrix2cd U1 = this->U1(s, k);
    std::complex<double> zTmp = phi0k.adjoint() * U1 * phi0k;
    return zTmp;

}

std::vector<std::pair<int, std::complex<double>>> producer::kAndSortZ(const double&s) {
    std::vector<std::complex<double>> zAllAtQ;
    for (int k = 0; k < this->CON.N + 1; k++) {
        //double dk=2*M_PI/(double)this->CON.N;
        zAllAtQ.push_back(this->calZ(s, (double) k*this->CON.dk));
    }

    std::vector<int> order(this->CON.N);
    int k0 = 0;
    std::iota(order.begin(), order.end(), k0++);
    std::sort(order.begin(), order.end(), [&](int i, int j) { return std::abs(zAllAtQ[i]) <= std::abs(zAllAtQ[j]); });
    std::vector<std::pair<int, std::complex<double>>> rst;
    for (const auto &elem:order) {
        rst.emplace_back(std::make_pair(elem, zAllAtQ[elem]));
    }
    return rst;
}

void producer::print123(const double&s) {
    std::vector<std::pair<int,std::complex<double>>> rst=this->kAndSortZ(s);
    std::cout<<"k = "<<(double)rst[0].first*this->CON.dk<<", z = "<<rst[0].second<<" ,abs(z) = "<<std::abs(rst[0].second)<<std::endl;
    std::cout<<"k = "<<(double)rst[1].first*this->CON.dk<<", z = "<<rst[1].second<<" ,abs(z) = "<<std::abs(rst[1].second)<<std::endl;
    std::cout<<"k = "<<(double)rst[2].first*this->CON.dk<<", z = "<<rst[2].second<<" ,abs(z) = "<<std::abs(rst[2].second)<<std::endl;



}

