#include"potentials.hpp"
#include"gk.hpp"
#include<iostream>

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

/*  Calculates (Z - Zt). No squaring here. Handle that inside whatever function 
    calls cardyError upstream. */
mpf_class cardyError(float ccharge, float beta, boost::multi_array<mpf_class, 1>& evals){
    mpfr_t ccharger, betar, betapr, zero;
    mpfr_init(ccharger);
    mpfr_set_flt(ccharger, ccharge, MPFR_RNDN);
    mpfr_init(betar);
    mpfr_set_flt(betar, beta, MPFR_RNDN);
    mpfr_init(betapr);
    mpfr_set_flt(betapr, 4*M_PI*M_PI, MPFR_RNDN);
    mpfr_div(betapr, betapr, betar, MPFR_RNDN);
    mpfr_init(zero);
    mpfr_set_str(zero, "0", 10, MPFR_RNDN);
    
    // This implements Z = sum(exp(- lambda beta)) and Zt = sum(exp(- lambda beta'))
    mpfr_t Z;
    mpfr_init_set_str(Z, "1", 10, MPFR_RNDN); // Initialise Z to 1 in base 10
    mpfr_t temp;
    mpfr_init(temp);

    mpfr_t Zt;
    mpfr_init_set_str(Zt, "1", 10, MPFR_RNDN); // Initialise Zt to 1 in base 10
    mpfr_t tempt;
    mpfr_init(tempt);  

    mpf_t debug;
    mpf_init(debug);

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
    mpfr_set_flt(vac, ccharge/(12), MPFR_RNDN);
    mpfr_mul(vac, vac, betar, MPFR_RNDN);
    mpfr_exp(vac,vac, MPFR_RNDN);

    mpfr_t vacCrossed;
    mpfr_init(vacCrossed);
    mpfr_set_flt(vacCrossed, ccharge/(12), MPFR_RNDN);
    mpfr_mul(vacCrossed, vacCrossed, betapr, MPFR_RNDN);
    mpfr_exp(vacCrossed, vacCrossed, MPFR_RNDN);

    mpfr_mul(Z, Z, vac, MPFR_RNDN);
    mpfr_mul(Zt, Zt, vacCrossed, MPFR_RNDN);

    // This implements Z - Zt
    mpfr_t error;
    mpfr_init(error);
    mpfr_set(error, Z, MPFR_RNDN);
    mpfr_sub(error, error, Zt, MPFR_RNDN);

    mpf_t returnValue;
    mpf_init(returnValue);
    mpfr_get_f(returnValue, error, MPFR_RNDN);
    mpf_class returnValueC(returnValue);

    mpfr_clear(ccharger);
    mpfr_clear(betar);
    mpfr_clear(betapr);
    mpfr_clear(zero);
    mpfr_clear(Z);
    mpfr_clear(Zt);
    mpfr_clear(error);
    mpfr_clear(vac);
    mpfr_clear(vacCrossed);
    mpf_clear(returnValue);

    return returnValueC;
}

/* Overloaded version of cardyError with mpfr_t input for beta */
mpf_class cardyError(float ccharge, mpfr_t beta, boost::multi_array<mpf_class, 1>& evals){
    mpfr_t ccharger, betar, betapr, zero;
    mpfr_init(ccharger);
    mpfr_set_flt(ccharger, ccharge, MPFR_RNDN);
    mpfr_init(betar);
    mpfr_set(betar, beta, MPFR_RNDN);
    mpfr_init(betapr);
    mpfr_set_flt(betapr, 4*M_PI*M_PI, MPFR_RNDN);
    mpfr_div(betapr, betapr, betar, MPFR_RNDN);
    mpfr_init(zero);
    mpfr_set_str(zero, "0", 10, MPFR_RNDN);
    
    // This implements Z = sum(exp(- lambda beta)) and Zt = sum(exp(- lambda beta'))
    mpfr_t Z;
    mpfr_init_set_str(Z, "1", 10, MPFR_RNDN); // Initialise Z to 1 in base 10
    mpfr_t temp;
    mpfr_init(temp);

    mpfr_t Zt;
    mpfr_init_set_str(Zt, "1", 10, MPFR_RNDN); // Initialise Zt to 1 in base 10
    mpfr_t tempt;
    mpfr_init(tempt);  

    mpf_t debug;
    mpf_init(debug);

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
    mpfr_set_flt(vac, ccharge/(12), MPFR_RNDN);
    mpfr_mul(vac, vac, betar, MPFR_RNDN);
    mpfr_exp(vac,vac, MPFR_RNDN);

    mpfr_t vacCrossed;
    mpfr_init(vacCrossed);
    mpfr_set_flt(vacCrossed, ccharge/(12), MPFR_RNDN);
    mpfr_mul(vacCrossed, vacCrossed, betapr, MPFR_RNDN);
    mpfr_exp(vacCrossed, vacCrossed, MPFR_RNDN);

    mpfr_mul(Z, Z, vac, MPFR_RNDN);
    mpfr_mul(Zt, Zt, vacCrossed, MPFR_RNDN);

    // This implements Z - Zt
    mpfr_t error;
    mpfr_init(error);
    mpfr_set(error, Z, MPFR_RNDN);
    mpfr_sub(error, error, Zt, MPFR_RNDN);

    mpf_t returnValue;
    mpf_init(returnValue);
    mpfr_get_f(returnValue, error, MPFR_RNDN);
    mpf_class returnValueC(returnValue);

    mpfr_clear(ccharger);
    mpfr_clear(betar);
    mpfr_clear(betapr);
    mpfr_clear(zero);
    mpfr_clear(Z);
    mpfr_clear(Zt);
    mpfr_clear(error);
    mpfr_clear(vac);
    mpfr_clear(vacCrossed);
    mpf_clear(returnValue);

    return returnValueC;
}

//  Returns exp(-(Z-Zt)**2) by calling cardyError for Z-Zt
mpf_class expV(float ccharge, float beta, int dim, boost::multi_array<mpf_class, 1>& evals){
    mpfr_t dimN;
    mpfr_init_set_ui(dimN, dim, MPFR_RNDN);

    mpf_class error = cardyError(ccharge, beta, evals);

    mpfr_t potential, pdf; 
    mpfr_init_set_str(pdf, "0", 10, MPFR_RNDN); // set pdf to 0 in base 10
    mpfr_init_set_f(potential, error.get_mpf_t(), MPFR_RNDN); // set pot = Z - Zt
    mpfr_add(potential, potential, dimN, MPFR_RNDN); // pot = Z - Zt + N
    mpfr_mul(potential, potential, potential, MPFR_RNDN); // pot = (Z - Zt + N) ** 2
    mpfr_sub(pdf, pdf, potential, MPFR_RNDN); // sets pdf = - (Z - Zt + N)**2;
    mpfr_exp(pdf, pdf, MPFR_RNDN); // exp(- (Z - Zt + N)**2) 
    // mpfr_printf ("potential is %.60Rf\n", potential);


    mpf_t returnValue;
    mpf_init(returnValue);
    mpfr_get_f(returnValue, pdf, MPFR_RNDN);
    
    // TODO - Is there a better way to achieve the following? For some reason
    // mpfr_clears() does not work and gives an error of undefined function
    mpfr_clear(potential);
    mpfr_clear(dimN);
    mpfr_clear(pdf);
    return mpf_class(returnValue);

    // No need to call deconstructor of error. Automatically happens when 
    // exiting scope
}


/*  Integrates (Z - Zt + N)**2 in a range between betamin betamax using the 
    trapezoid method and returns exp(-integral)
*/
mpf_class expV_int_trap(float ccharge, float betamin, float betamax, unsigned int dim, boost::multi_array<mpf_class, 1>& evals, int num_steps) {

    mpfr_t dimN;
    mpfr_init_set_ui(dimN, dim, MPFR_RNDN);

    mpf_class error;
    mpf_t step_contrib;
    mpf_init(step_contrib);

    mpf_class integrated_result = 0;
    mpfr_t potential, pdf; 
    mpfr_init(potential);
    mpfr_init(pdf);

    float step = (betamax - betamin) / num_steps;
    for (int i = 0; i <= num_steps; ++i) {
        error = cardyError(ccharge, betamin + i * step, evals);

        mpfr_set_str(pdf, "0", 10, MPFR_RNDN); // set pdf to 0 in base 10
        mpfr_set_f(potential, error.get_mpf_t(), MPFR_RNDN); // set pot = Z - Zt
        // mpfr_add(potential, potential, dimN, MPFR_RNDN); // pot = Z - Zt + N
        mpfr_mul(potential, potential, potential, MPFR_RNDN); // pot = (Z - Zt + N) ** 2

        mpfr_get_f(step_contrib, potential, MPFR_RNDN);

        if (i == 0 || i == num_steps) {
            integrated_result += mpf_class(step_contrib) * 0.5 * step;  // Edge cases (0.5 weight for first/last step)
        } 
        else {
            integrated_result += mpf_class(step_contrib) * step;  // Middle steps (full weight)
        }
    }
    mpfr_t integral_result;
    mpfr_init_set_f(integral_result, integrated_result.get_mpf_t(), MPFR_RNDN);
    mpfr_set_str(pdf, "0", 10, MPFR_RNDN);
    mpfr_sub(pdf, pdf, integral_result, MPFR_RNDN); // sets pdf = - integral
    mpfr_exp(pdf, pdf, MPFR_RNDN); // exp(- integral) 


    // Convert the result back to mpf_class
    mpf_t returnValue;
    mpf_init(returnValue);
    mpfr_get_f(returnValue, pdf, MPFR_RNDN);


    return mpf_class(returnValue);
}

/*  Integrates (Z - Zt + N)**2 with Hermite-Gauss quadrature and then returns
    exp(- N * integral)
*/
mpf_class expV_GK(float ccharge, int order, unsigned int dim, boost::multi_array<mpf_class, 1>& evals) {
    /* int order = The number of points where the integral is evaluated. Also 
                    equal to the number of zeros of the Hermite polynomial
    
    */

    // Variables necessary to hold the answers
    mpfr_t integral, mpfr_xi, mpfr_wi, integrand;
    mpfr_init_set_str(integral,"0", 10, MPFR_RNDN);
    mpfr_init(mpfr_xi);
    mpfr_init(mpfr_wi);
    mpfr_init(integrand);
    mpf_class error;
    mpfr_t dimN;
    mpfr_init_set_ui(dimN, dim, MPFR_RNDN);
    mpfr_t zero;
    mpfr_init_set_str(zero, "0", 10, MPFR_RNDN);

    #ifdef DEBUG
    mpf_t debug;
    mpf_init(debug);
    #endif

    // Code to fetch the abcissa and weights
    const auto& xi = gk_abcissa(order);
    const auto& wi = gk_weights(order);

    for (size_t i = 0; i < xi.size(); ++i) {
        const char* x = xi[i].c_str();
        const char* w = wi[i].c_str();

        if (x[0] == '-'){
            continue;
        }

        mpfr_set_str(mpfr_xi, x, 10, MPFR_RNDN);
        mpfr_add_d(mpfr_xi, mpfr_xi, 2*M_PI, MPFR_RNDN); // Offset by 2pi since the integral is from 2pi to infinity
        mpfr_set_str(mpfr_wi, w, 10, MPFR_RNDN);

        // Code to calculate the potential at each abcissa
        error = cardyError(ccharge, mpfr_xi, evals);

        #ifdef DEBUG
        std::cout<<error<<std::endl;
        #endif

        mpfr_set_f(integrand, error.get_mpf_t(), MPFR_RNDN); // set integrand = Z - Zt at beta = xi
        // mpfr_add(integrand, integrand, dimN, MPFR_RNDN); // integrand = Z - Zt + N
        mpfr_mul(integrand, integrand, integrand, MPFR_RNDN); // integrand = (Z - Zt + N) ** 2

        mpfr_mul(integrand, integrand, mpfr_wi, MPFR_RNDN); // weight the integrand
        mpfr_add(integral, integral, integrand, MPFR_RNDN); // add to integral
    }
    // Exponentiate the integrated result

    #ifdef DEBUG
    mpfr_get_f(debug, integral, MPFR_RNDN);
    std::cout<<mpf_class(debug)<<std::endl;
    #endif

    mpfr_sub(integral, zero, integral, MPFR_RNDN); // integral = -integral
    mpfr_mul(integral, integral, dimN, MPFR_RNDN); // integral = - N * integral
    mpfr_exp(integral, integral, MPFR_RNDN); // integral = exp(- N * integral)

    // Return
    mpf_t returnValue;
    mpf_init(returnValue);
    mpfr_get_f(returnValue, integral, MPFR_RNDN);

    mpfr_clear(integral);
    mpfr_clear(mpfr_xi);
    mpfr_clear(mpfr_wi);
    mpfr_clear(integrand);
    mpfr_clear(dimN);
    mpfr_clear(zero);

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
    mpfr_get_f(returnValue, Z, MPFR_RNDN);  // Convert mpfr_t Z to mpf_t

    // Clear memory used by mpfr_t variables
    mpfr_clear(temp);
    mpfr_clear(Z);
    mpfr_clear(betaRegr);
    mpfr_clear(Ewallr);

    mpf_class returnValueC(returnValue);
    mpf_clear(returnValue);

    // Return the final result as mpf_class
    return returnValueC;
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

mpf_class Cardy_int_trap(float ccharge, float betamin, float betamax, int dim, boost::multi_array<mpf_class, 1> evals){
    return vanderMonde(evals) * expV_int_trap(ccharge, betamin, betamax, dim, evals, 100)* gaussianwall(1, evals, 200); // 100 is a placeholder
}

mpf_class Cardy_int_GK(float ccharge, int order, int dim, float Ewall, boost::multi_array<mpf_class, 1> evals){
    return vanderMonde(evals) * expV_GK(ccharge, order, dim, evals) * gaussianwall(1, evals, Ewall);//; // 100 is a placeholder
}