#include"boost/multi_array.hpp"
#include<gmp.h>
#include<gmpxx.h>
#include<mpfr.h>

struct cardyParams {
    float ccharge;
    float beta;
};

mpf_class vanderMonde(boost::multi_array<mpf_class, 1> evals){
    int len = evals.size();
    mpf_class vm = 1;
    for (int i = 0; i < len; i++){
        for (int j = i + 1; j< len; j++){
            vm *= (evals[i] - evals[j]);
        }
    }
    vm *= vm;
    return vm;
}

mpf_class expV(float ccharge, float beta, boost::multi_array<mpf_class, 1> evals){
    mpfr_t ccharger, betar, betapr, zero;
    mpfr_init(ccharger);
    mpfr_set_flt(ccharger, ccharge, MPFR_RNDN);
    mpfr_init(betar);
    mpfr_set_flt(betar, beta, MPFR_RNDN);
    mpfr_init(betapr);
    mpfr_set_flt(betapr, 4*M_PI*M_PI/beta, MPFR_RNDN);
    mpfr_init(zero);
    mpfr_set_str(zero, "0", 10, MPFR_RNDN);
    
    // This implements Z = sum(exp(lambda beta)) and Zt = sum(exp(lambda beta,))
    mpfr_t Z;
    mpfr_init_set_str(Z, "1", 10, MPFR_RNDN); // Initialise Z to 1 in base 10
    mpfr_t temp;
    mpfr_init(temp);

    mpfr_t Zt;
    mpfr_init_set_str(Zt, "1", 10, MPFR_RNDN); // Initialise Zt to 1 in base 10
    mpfr_t tempt;
    mpfr_init(tempt);

    mpfr_t dim; // TODO - CHECK that N appears in front of Z-Zt
    mpfr_init_set_str(dim, "100", 10, MPFR_RNDN);    

    for (auto x: evals){
        mpfr_set_f(temp, x.get_mpf_t(), MPFR_RNDN);
        mpfr_sub(temp, zero, temp, MPFR_RNDN);
        mpfr_mul(temp, temp, betar, MPFR_RNDN);
        mpfr_exp(temp, temp, MPFR_RNDN);
        mpfr_add(Z, Z, temp, MPFR_RNDN);

        mpfr_set_f(tempt, x.get_mpf_t(), MPFR_RNDN);
        mpfr_sub(tempt, zero, tempt, MPFR_RNDN);
        mpfr_mul(tempt, tempt, betapr, MPFR_RNDN);
        mpfr_exp(tempt, tempt, MPFR_RNDN);
        mpfr_add(Zt, Zt, tempt, MPFR_RNDN);
    }
    // mpfr_printf ("Z is %.60Rf\n", Z);
    // mpfr_printf ("Zt is %.60Rf\n", Zt);
    // mpfr_clears(temp,tempt);

    mpfr_t vac;
    mpfr_init(vac);
    mpfr_set_flt(vac, ccharge*beta/(12), MPFR_RNDN);
    mpfr_t vacCrossed;
    mpfr_init(vacCrossed);
    mpfr_set_flt(vacCrossed, ccharge/(12) * 4*M_PI*M_PI/beta, MPFR_RNDN);

    mpfr_mul(Z, Z, vac, MPFR_RNDN);
    mpfr_mul(Zt, Zt, vacCrossed, MPFR_RNDN);
    
    

    // This implements exp((Z - Zt)**2)
    mpfr_t potential;
    mpfr_init(potential);
    mpfr_set(potential, Z, MPFR_RNDN);
    mpfr_sub(potential, potential, Zt, MPFR_RNDN);
    mpfr_add(potential, potential, dim, MPFR_RNDN); // Z - Zt + N
    mpfr_mul(potential, potential, potential, MPFR_RNDN); // (Z - Zt + N)**2
    mpfr_exp(potential, potential, MPFR_RNDN); // exp((Z - Zt + N)**2) 
    // mpfr_printf ("potential is %.60Rf\n", potential);


    mpf_t returnValue;
    mpf_init(returnValue);
    mpfr_get_f(returnValue, potential, MPFR_RNDN);
    
    // TODO - Is there a better way to achieve the following? For some reason
    // mpfr_clears() does not work and gives an error of undefined function
    mpfr_clear(ccharger);
    mpfr_clear(betar);
    mpfr_clear(betapr);
    mpfr_clear(zero);
    mpfr_clear(Z);
    mpfr_clear(Zt);
    mpfr_clear(potential);
    return mpf_class(returnValue);
}

mpf_class gaussian(float betaReg, boost::multi_array<mpf_class, 1> evals){
    // This implements Z = exp(-betaReg * sum(x**2))
    mpfr_t betaRegr;
    mpfr_init(betaRegr);
    mpfr_set_flt(betaRegr, betaReg, MPFR_RNDN);

    mpfr_t Z;
    mpfr_init_set_str(Z, "0", 10, MPFR_RNDN); // Initialise Z to 0 in base 10
    
    mpfr_t temp;
    mpfr_init(temp);
    for (auto x: evals){
        mpfr_set_f(temp, x.get_mpf_t(), MPFR_RNDN); // eval[i]
        mpfr_mul(temp, temp, temp, MPFR_RNDN); // eval[i]**2
        mpfr_sub(Z, Z, temp, MPFR_RNDN); // Z - eval[i]**2
    }
    mpfr_clear(temp);
    mpfr_mul(Z, betaRegr, Z, MPFR_RNDN);
    mpfr_exp(Z, Z, MPFR_RNDN); // Z <- exp(Z)
    
    mpf_t returnValue;
    mpf_init(returnValue);
    mpfr_get_f(returnValue, Z, MPFR_RNDN); // Need to change to potential
    mpfr_clear(Z);
    return mpf_class(returnValue);
}

mpf_class wigner(float betaReg, boost::multi_array<mpf_class, 1> evals){
    return vanderMonde(evals) * gaussian(betaReg, evals);
}

mpf_class cardy(float ccharge, float beta, boost::multi_array<mpf_class, 1> evals){
    return vanderMonde(evals) * expV(ccharge, beta, evals);
}

mpf_class dampedCardy(float ccharge, float beta, float betaReg, boost::multi_array<mpf_class, 1> evals){
    return vanderMonde(evals) * expV(ccharge, beta, evals);
}