#include"boost/multi_array.hpp"
#include<gmp.h>
#include<gmpxx.h>
#include<mpfr.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <vector>
#include <cmath>

using namespace boost::multiprecision;
using namespace boost::math::quadrature;

struct cardyParams {
    float ccharge;
    float beta;
};

mpf_class vanderMonde(boost::multi_array<mpf_class, 1>& evals){
    int len = evals.size();
    mpf_class vm = 1;
    for (int i = 0; i < len; i++){
        for (int j = i + 1; j< len; j++){
            vm *= (evals[i] - evals[j]);
        }
        vm *= (evals[i]);
    }
    vm *= vm;
    return vm;
}

mpf_class expV(float ccharge, float beta, int dim, boost::multi_array<mpf_class, 1>& evals){
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

    mpfr_t dimN; // TODO - CHECK that N appears in front of Z-Zt
    mpfr_init_set_ui(dimN, dim, MPFR_RNDN);    

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
    mpfr_t potential, pdf;
    mpfr_init(potential);
    mpfr_init_set_str(pdf, "0", 10, MPFR_RNDN);
    mpfr_set(potential, Z, MPFR_RNDN);
    mpfr_sub(potential, potential, Zt, MPFR_RNDN);
    //mpfr_add(potential, potential, dimN, MPFR_RNDN); // Z - Zt + N
    mpfr_mul(potential, potential, potential, MPFR_RNDN); // (Z - Zt + N)**2
    mpfr_sub(pdf, pdf, potential, MPFR_RNDN); // sets pdf = - (Z - Zt + N)**2;
    mpfr_exp(potential, potential, MPFR_RNDN); // exp(- (Z - Zt + N)**2) 
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
    mpfr_clear(dimN);
    mpfr_clear(pdf);
    return mpf_class(returnValue);
}


//integrates (Z-Zt)**2 on a beta range between betamin betamax using the trapezoid method
mpf_class expV_int_trap(float ccharge, float betamin, float betamax, int dim, boost::multi_array<mpf_class, 1>& evals, int num_steps) {
    // Initialize variables as before
    mpfr_t ccharger, zero;
    mpfr_init(ccharger);
    mpfr_set_flt(ccharger, ccharge, MPFR_RNDN);

    mpfr_init(zero);
    mpfr_set_str(zero, "0", 10, MPFR_RNDN);

    mpfr_t Z, Zt, temp, tempt, vac, vacCrossed, dimN, betar, betapr, potential, pdf;
    mpfr_init_set_str(Z, "1", 10, MPFR_RNDN); // Initialize Z to 1
    mpfr_init_set_str(Zt, "1", 10, MPFR_RNDN); // Initialize Zt to 1
    mpfr_init(temp);
    mpfr_init(tempt);
    mpfr_init(vac);
    mpfr_init(vacCrossed);
    mpfr_init_set_ui(dimN, dim, MPFR_RNDN);  // Initialize dimN to `dim`
    mpfr_init(potential);
    mpfr_init(pdf);
    mpfr_init(betar);
    mpfr_init(betapr);

    // Initialize the return variable for the accumulated integral
    mpf_class integrated_result = 0;

    // Numerical integration over beta using trapezoidal rule
    float step = (betamax - betamin) / num_steps;
    for (int i = 0; i <= num_steps; ++i) {
        float beta = betamin + i * step;

        // Set `betar` and `betapr` as the current value of beta and 4π²/beta
        mpfr_set_flt(betar, beta, MPFR_RNDN);
        mpfr_set_flt(betapr, 4 * M_PI * M_PI / beta, MPFR_RNDN);

        // Reset Z and Zt for each beta value
        mpfr_set_str(Z, "1", 10, MPFR_RNDN);
        mpfr_set_str(Zt, "1", 10, MPFR_RNDN);

        // Compute Z and Zt as sum(exp(lambda * beta)) and sum(exp(lambda * 4π² / beta))
        for (auto x : evals) {
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

        // Compute vac and vacCrossed for the current beta
        mpfr_set_flt(vac, ccharge * beta / 12, MPFR_RNDN);
        mpfr_set_flt(vacCrossed, ccharge / 12 * 4 * M_PI * M_PI / beta, MPFR_RNDN);

        mpfr_mul(Z, Z, vac, MPFR_RNDN);
        mpfr_mul(Zt, Zt, vacCrossed, MPFR_RNDN);

        // Compute (Z - Zt + N)^2
        mpfr_set(potential, Z, MPFR_RNDN);
        mpfr_sub(potential, potential, Zt, MPFR_RNDN);
        mpfr_add(potential, potential, dimN, MPFR_RNDN);  // (Z - Zt + N)
        mpfr_mul(potential, potential, potential, MPFR_RNDN);  // (Z - Zt + N)^2

        // Compute - (Z - Zt + N)^2 and exp of the negative term
        mpfr_sub(pdf, zero, potential, MPFR_RNDN);  // - (Z - Zt + N)^2

        // Accumulate the potential for numerical integration (trapezoidal rule)
        mpf_t step_contrib;
        mpf_init(step_contrib);
        mpfr_get_f(step_contrib, pdf, MPFR_RNDN);  // Get the current potential contribution

        if (i == 0 || i == num_steps) {
            integrated_result += mpf_class(step_contrib) * 0.5 * step;  // Edge cases (0.5 weight for first/last step)
        } else {
            integrated_result += mpf_class(step_contrib) * step;  // Middle steps (full weight)
        }

        // Clear local variables for the loop
        mpf_clear(step_contrib);
    }

    // Exponentiate the integrated result: exp( integrated_result )
    mpfr_t final_exp_result;
    mpfr_init(final_exp_result);
    mpfr_set_f(final_exp_result, integrated_result.get_mpf_t(), MPFR_RNDN);  // Set the final integral
    mpfr_exp(final_exp_result, final_exp_result, MPFR_RNDN);  // exp(integrated_result)

    // Convert the result back to mpf_class
    mpf_t returnValue;
    mpf_init(returnValue);
    mpfr_get_f(returnValue, final_exp_result, MPFR_RNDN);

    // Clear mpfr variables
    mpfr_clear(ccharger);
    mpfr_clear(zero);
    mpfr_clear(Z);
    mpfr_clear(Zt);
    mpfr_clear(temp);
    mpfr_clear(tempt);
    mpfr_clear(vac);
    mpfr_clear(vacCrossed);
    mpfr_clear(dimN);
    mpfr_clear(potential);
    mpfr_clear(pdf);
    mpfr_clear(betar);
    mpfr_clear(betapr);
    mpfr_clear(final_exp_result);

    return mpf_class(returnValue);
}

//implements the exp of integral (Z-Zt)**2 using the gaussian-kronrod quadratures, with maximum subdivision steps = maxdepth 
mpf_class expV_int_gk(float ccharge, float betamin, float betamax, int dim, 
                      boost::multi_array<mpf_class, 1>& evals, int max_depth) {
    // Initialize variables
    mpfr_t ccharger, zero;
    mpfr_init(ccharger);
    mpfr_set_flt(ccharger, ccharge, MPFR_RNDN);

    mpfr_init(zero);
    mpfr_set_str(zero, "0", 10, MPFR_RNDN);

    mpfr_t Z, Zt, temp, tempt, vac, vacCrossed, dimN, betar, betapr, potential;
    mpfr_init_set_str(Z, "1", 10, MPFR_RNDN); // Initialize Z to 1
    mpfr_init_set_str(Zt, "1", 10, MPFR_RNDN); // Initialize Zt to 1
    mpfr_init(temp);
    mpfr_init(tempt);
    mpfr_init(vac);
    mpfr_init(vacCrossed);
    mpfr_init_set_ui(dimN, dim, MPFR_RNDN);  // Initialize dimN to dim
    mpfr_init(potential);
    mpfr_init(betar);
    mpfr_init(betapr);

    // Function for Gauss-Kronrod integration
    auto integrand = [&](double beta) {
        // Set betar and betapr as the current value of beta and 4π²/beta
        mpfr_set_flt(betar, beta, MPFR_RNDN);
        mpfr_set_flt(betapr, 4 * M_PI * M_PI / beta, MPFR_RNDN);

        // Reset Z and Zt for each beta value
        mpfr_set_str(Z, "1", 10, MPFR_RNDN);
        mpfr_set_str(Zt, "1", 10, MPFR_RNDN);

        // Compute Z and Zt as sum(exp(lambda * beta)) and sum(exp(lambda * 4π² / beta))
        for (auto x : evals) {
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

        // Compute vac and vacCrossed for the current beta
        mpfr_set_flt(vac, ccharge * beta / 12, MPFR_RNDN);
        mpfr_set_flt(vacCrossed, ccharge / 12 * 4 * M_PI * M_PI / beta, MPFR_RNDN);

        mpfr_mul(Z, Z, vac, MPFR_RNDN);
        mpfr_mul(Zt, Zt, vacCrossed, MPFR_RNDN);

        // Compute (Z - Zt + N)^2
        mpfr_set(potential, Z, MPFR_RNDN);
        mpfr_sub(potential, potential, Zt, MPFR_RNDN);
        //for now avoid this regularization, needed to have finite int when betamax=infty
        //mpfr_add(potential, potential, dimN, MPFR_RNDN);  // (Z - Zt + N)
        mpfr_mul(potential, potential, potential, MPFR_RNDN);  // (Z - Zt + N)^2

        // Get the current potential contribution
        double result = mpfr_get_d(potential, MPFR_RNDN);

        return result;
    };

    // Perform Gauss-Kronrod integration over [betamin, betamax]
    double integral_result = gauss_kronrod<double, 15>::integrate(integrand, betamin, betamax, max_depth);

    // Exponentiate the integrated result: exp( integral_result )
    mpfr_t final_exp_result;
    mpfr_init(final_exp_result);
    mpfr_set_d(final_exp_result, integral_result, MPFR_RNDN);
    mpfr_exp(final_exp_result, final_exp_result, MPFR_RNDN);  // exp(integral_result)

    // Convert the result back to mpf_class
    mpf_t returnValue;
    mpf_init(returnValue);
    mpfr_get_f(returnValue, final_exp_result, MPFR_RNDN);

    // Clear mpfr variables
    mpfr_clear(ccharger);
    mpfr_clear(zero);
    mpfr_clear(Z);
    mpfr_clear(Zt);
    mpfr_clear(temp);
    mpfr_clear(tempt);
    mpfr_clear(vac);
    mpfr_clear(vacCrossed);
    mpfr_clear(dimN);
    mpfr_clear(potential);
    mpfr_clear(betar);
    mpfr_clear(betapr);
    mpfr_clear(final_exp_result);

    return mpf_class(returnValue);
}




mpf_class expVHigh(float ccharge, float beta, int dim, boost::multi_array<mpf_class, 1>& evals){
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

    mpfr_t dimN; // TODO - CHECK that N appears in front of Z-Zt
    mpfr_init_set_ui(dimN, dim, MPFR_RNDN);    

    for (auto x: evals){
        mpfr_set_f(temp, x.get_mpf_t(), MPFR_RNDN);
        mpfr_sub(temp, zero, temp, MPFR_RNDN);
        mpfr_mul(temp, temp, betar, MPFR_RNDN);
        mpfr_exp(temp, temp, MPFR_RNDN);
        mpfr_add(Z, Z, temp, MPFR_RNDN);

        //mpfr_set_f(tempt, x.get_mpf_t(), MPFR_RNDN);
        //mpfr_sub(tempt, zero, tempt, MPFR_RNDN);
        //mpfr_mul(tempt, tempt, betapr, MPFR_RNDN);
        //mpfr_exp(tempt, tempt, MPFR_RNDN);
        //mpfr_add(Zt, Zt, tempt, MPFR_RNDN);
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
    mpfr_t potential, pdf;
    mpfr_init(potential);
    mpfr_init_set_str(pdf, "0", 10, MPFR_RNDN);
    mpfr_set(potential, Z, MPFR_RNDN);
    mpfr_sub(potential, potential, Zt, MPFR_RNDN);
    mpfr_add(potential, potential, dimN, MPFR_RNDN); // Z - Zt + N
    mpfr_mul(potential, potential, potential, MPFR_RNDN); // (Z - Zt + N)**2
    mpfr_sub(pdf, pdf, potential, MPFR_RNDN); // sets pdf = - (Z - Zt + N)**2;
    mpfr_exp(potential, potential, MPFR_RNDN); // exp(- (Z - Zt + N)**2) 
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
    mpfr_clear(dimN);
    mpfr_clear(pdf);
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





mpf_class gaussianwall(float betaReg, const boost::multi_array<mpf_class, 1>& evals, float Ewall) {
    // Initialize mpfr_t to hold betaReg with arbitrary precision
    mpfr_t betaRegr;
    mpfr_init(betaRegr);
    mpfr_set_flt(betaRegr, betaReg, MPFR_RNDN);  // Set betaRegr from float betaReg

    // Initialize mpfr_t to hold Ewall with arbitrary precision
    mpfr_t Ewallr;
    mpfr_init(Ewallr);
    mpfr_set_flt(Ewallr, Ewall, MPFR_RNDN);  // Set Ewallr from float Ewall

    // Initialize Z to 0 (the result accumulator)
    mpfr_t Z;
    mpfr_init_set_str(Z, "0", 10, MPFR_RNDN);  // Initialise Z to 0 in base 10

    // Temporary variable for calculations
    mpfr_t temp;
    mpfr_init(temp);

    // Loop through each eigenvalue in the evals array
    for (size_t i = 0; i < evals.size(); ++i) {
        // Get the eigenvalue `x`
        mpfr_set_f(temp, evals[i].get_mpf_t(), MPFR_RNDN);  // Set temp to the current eigenvalue

        // Check if x > Ewall
        if (mpfr_cmp(temp, Ewallr) > 0) {
            // x > Ewall, so we apply the Gaussian calculation
            
            // Subtract Ewall from the eigenvalue (temp = x - Ewall)
            mpfr_sub(temp, temp, Ewallr, MPFR_RNDN);
            
            // Square the result (temp = (x - Ewall)^2)
            mpfr_mul(temp, temp, temp, MPFR_RNDN);

            // Subtract temp from Z (Z = Z - (x - Ewall)^2)
            mpfr_sub(Z, Z, temp, MPFR_RNDN);
        }
    }

    // Multiply Z by betaReg (Z = betaReg * Z)
    mpfr_mul(Z, betaRegr, Z, MPFR_RNDN);

    // Exponentiate Z (Z = exp(Z))
    mpfr_exp(Z, Z, MPFR_RNDN);

    // Convert Z to mpf_class to return the result
    mpf_t returnValue;
    mpf_init(returnValue);
    mpfr_get_f(returnValue, Z, MPFR_RNDN);  // Convert mpfr_t Z to mpf_class

    // Clear memory used by mpfr_t variables
    mpfr_clear(temp);
    mpfr_clear(Z);
    mpfr_clear(betaRegr);
    mpfr_clear(Ewallr);

    // Return the final result as mpf_class
    return mpf_class(returnValue);
}


mpf_class wigner(float betaReg, boost::multi_array<mpf_class, 1> evals){
    return vanderMonde(evals) * gaussian(betaReg, evals);
}

mpf_class cardy(float ccharge, float beta, int dim, boost::multi_array<mpf_class, 1>& evals){
    return vanderMonde(evals) * expV(ccharge, beta, dim, evals);
}

mpf_class dampedCardy(float ccharge, float beta, float betaReg, boost::multi_array<mpf_class, 1> evals){
    return vanderMonde(evals) * expV(ccharge, beta, 100, evals); // 100 is a placeholder
}
mpf_class CardyGwall(float ccharge, float beta, int dim, boost::multi_array<mpf_class, 1> evals){
    return vanderMonde(evals) * expV(ccharge, beta, dim, evals)* gaussianwall(1, evals, 200); // 100 is a placeholder
}

mpf_class CardyHigh(float ccharge, float beta, int dim, boost::multi_array<mpf_class, 1> evals){
    return vanderMonde(evals) * expVHigh(ccharge, beta, dim, evals)* gaussianwall(1, evals, 200); // 100 is a placeholder
}

mpf_class Cardy_int_trap(float ccharge, float betamin, float betamax, int dim, boost::multi_array<mpf_class, 1> evals){
    return vanderMonde(evals) * expV_int_trap(ccharge, betamin, betamax, dim, evals,100)* gaussianwall(1, evals, 200); // 100 is a placeholder
}

mpf_class Cardy_int_gk(float ccharge, float betamin, float betamax, int dim, boost::multi_array<mpf_class, 1> evals){
    return vanderMonde(evals) * expV_int_gk(ccharge, betamin, betamax, dim, evals, 15)* gaussianwall(1, evals, 200); // 100 is a placeholder
}
