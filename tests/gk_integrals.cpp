#include "../src/potentials/gk.hpp"
#include "../src/potentials/potentials.hpp"
#include<gmp.h>
#include<gmpxx.h>
#include<mpfr.h>
#include<iostream>

#define PREC 500

#ifndef DEBUG
#define DEBUG 1
#endif

mpf_class gaussian_int(int order){
    const auto& xi = gk_abcissa(order);
    const auto& wi = gk_weights(order);

    mpfr_t integral, mpfr_xi, mpfr_wi, integrand;
    mpfr_init_set_str(integral,"0", 10, MPFR_RNDN);
    mpfr_init(mpfr_xi);
    mpfr_init(mpfr_wi);
    mpfr_init(integrand);
    for (size_t i = 0; i < xi.size(); ++i) {
        const char* x = xi[i].c_str();
        const char* w = wi[i].c_str();

        mpfr_set_str(mpfr_xi, x, 10, MPFR_RNDN);
        mpfr_set_str(mpfr_wi, w, 10, MPFR_RNDN);

        mpfr_set_str(integrand, "1", 10, MPFR_RNDN);

        mpfr_mul(integrand, integrand, mpfr_wi, MPFR_RNDN);
        mpfr_add(integral, integral, integrand, MPFR_RNDN);
    }
    mpf_t returnValue;
    mpf_init(returnValue);
    mpfr_get_f(returnValue, integral, MPFR_RNDN);
    mpfr_printf ("Integral is %.2000Rf\n", integral);

    mpfr_clear(integral);

    return mpf_class(returnValue);
}

int main(){
    mpf_set_default_prec(PREC);
    mpfr_set_default_prec(PREC);
    // std::cout<<gaussian_int(17)<<std::endl;

    boost::multi_array<mpf_class, 1> evals (boost::extents[40]);
    evals[0] = 0;
    evals[1] = 10;
    evals[2] = 20;
    evals[3] = 30;
    evals[4] = 40;
    evals[5] = 50;
    evals[6] = 60;
    evals[7] = 70;
    evals[8] = 80;
    evals[9] = 90;
    evals[10] = 100;
    evals[11] = 110;
    evals[12] = 120;
    evals[13] = 130;
    evals[14] = 140;
    evals[15] = 150;
    evals[16] = 160;
    evals[17] = 170;
    evals[18] = 180;
    evals[19] = 190;
    evals[20] = 200;
    evals[21] = 210;
    evals[22] = 220;
    evals[23] = 230;
    evals[24] = 240;
    evals[25] = 250;
    evals[26] = 260;
    evals[27] = 270;
    evals[28] = 280;
    evals[29] = 290;
    evals[30] = 300;
    evals[31] = 310;
    evals[32] = 320;
    evals[33] = 330;
    evals[34] = 340;
    evals[35] = 350;
    evals[36] = 360;
    evals[37] = 370;
    evals[38] = 380;
    evals[39] = 390;
    
    std::cout<<Cardy_int_GK(2, 18, 40, 370, evals)<<std::endl;

    return 0;
}