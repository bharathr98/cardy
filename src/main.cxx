#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<memory>
#include<gmp.h>
#include<gmpxx.h>
#include<mpfr.h>
#include<chrono>
#include"boost/multi_array.hpp"


#define PREC 500


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
    mpfr_t ccharger, betar, betapr;
    mpfr_init(ccharger);
    mpfr_set_flt(ccharger, ccharge, MPFR_RNDN);
    mpfr_init(betar);
    mpfr_set_flt(betar, beta, MPFR_RNDN);
    mpfr_init(betapr);
    mpfr_set_flt(betapr, 4*M_PI*M_PI/beta, MPFR_RNDN);
    
    // This implements Z = sum(exp(lambda beta)) and Zt = sum(exp(lambda beta,))
    mpfr_t Z;
    mpfr_init_set_str(Z, "0", 10, MPFR_RNDN); // Initialise Z to 0 in base 10
    mpfr_t temp;
    mpfr_init(temp);

    mpfr_t Zt;
    mpfr_init_set_str(Zt, "0", 10, MPFR_RNDN); // Initialise Z to 0 in base 10
    mpfr_t tempt;
    mpfr_init(tempt);

    for (auto x: evals){
        mpfr_set_f(temp, x.get_mpf_t(), MPFR_RNDN);
        mpfr_mul(temp, temp, betar, MPFR_RNDN);
        mpfr_exp(temp, temp, MPFR_RNDN);
        mpfr_add(Z, Z, temp, MPFR_RNDN);

        mpfr_set_f(tempt, x.get_mpf_t(), MPFR_RNDN);
        mpfr_mul(tempt, tempt, betapr, MPFR_RNDN);
        mpfr_exp(tempt, tempt, MPFR_RNDN);
        mpfr_add(Zt, Zt, tempt, MPFR_RNDN);
    }
    mpfr_printf ("Z is %.60Rf\n", Z);
    mpfr_printf ("Zt is %.60Rf\n", Zt);
    // mpfr_clears(temp,tempt);
    

    // This implements exp(Z - Zt)
    // mpfr_t potential;
    // mpfr_init(potential);
    // mpfr_set(potential, Z, MPFR_RNDN);
    // mpfr_sub(potential, potential, Zt, MPFR_RNDN);
    // mpfr_exp(potential, potential, MPFR_RNDN);
    // mpfr_printf ("potential is %.60Rf\n", potential);


    mpf_t returnValue;
    mpf_init(returnValue);
    mpfr_get_f(returnValue, Z, MPFR_RNDN); // Need to change to potential
    return mpf_class(returnValue);
}

mpf_class gaussian(boost::multi_array<mpf_class, 1> evals){
    // This implements Z = exp(-sum(x**2)) and Zt = sum(exp(lambda beta,))
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
    mpfr_exp(Z, Z, MPFR_RNDN); // Z <- exp(Z)
    
    mpf_t returnValue;
    mpf_init(returnValue);
    mpfr_get_f(returnValue, Z, MPFR_RNDN); // Need to change to potential
    mpfr_clear(Z);
    return mpf_class(returnValue);
}

mpf_class wigner(boost::multi_array<mpf_class, 1> evals){
    return vanderMonde(evals) * gaussian(evals);
}


unsigned long long int read_urandom()
{
	union {
		unsigned long long int value;
		char cs[sizeof(unsigned long long int)];
	} u;

	std::ifstream rfin("/dev/urandom");
	rfin.read(u.cs, sizeof(u.cs));
	rfin.close();

	return u.value;
}

int main(){
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    mpf_set_default_prec(PREC);
    mpfr_set_default_prec(PREC);
    std::cout<<"Precision is: "<<mpfr_get_default_prec()<<std::endl;
    gmp_randclass rand (gmp_randinit_default);

    mpf_class low, high, step_size;
    low = -5;
    high = 5;
    step_size = 0.1;

    // Create a 3D array that is maxTime x numWalkers x dim
    int maxTime = 100000;
    int numWalkers = 5;
    int dim = 50;
    typedef boost::multi_array<mpf_class, 3> array_type; // 3 here is the depth.
    boost::array<array_type::index, 3> shape = {{ maxTime, numWalkers, dim }}; // 3 here seems to do nothing as long as it is greater than depth. 
    array_type data(shape);
    
    // Randomly initialise the first time step
    for(int i = 0; i < numWalkers; i++){
        for(int j = 0; j < dim; j++){
            rand.seed(read_urandom());
            data[0][i][j] = (high - low)*rand.get_f() + low;
        }
    }

    // Evolve the walks
    typedef boost::multi_array<mpf_class, 1> step_type; // 1 here is the depth.
    boost::array<step_type::index, 1> stepShape = {{ dim }};
    step_type nextStep(stepShape);
    for(int time = 1; time < maxTime; time++){
        for(int i = 0; i < numWalkers; i++){
            for(int d = 0; d < dim; d++){
                nextStep[d] = data[time-1][i][d] + (2*step_size)*rand.get_f() - step_size;
            }
            mpf_class ratioPotential = wigner(nextStep)/wigner(data[time-1][i]);
            mpf_class decisionToss = rand.get_f();

            if (i == 0 && time % 500 == 0){
                std::cout<<"ratioPotential at time = "<<time<<" is = "<<ratioPotential<<std::endl;
            }

            if (decisionToss > ratioPotential){
                    data[time][i] = data[time - 1][i];
            }
            else{
                data[time][i] = nextStep;
            }
        }
    }

    if(dim == 1){
        for(int i = 0; i < numWalkers; i++){
            std::cout<<"ID "<<i;
            mpf_class sum = 0;
            for(int time = 0; time < maxTime; time++){
                sum += data[time][i][0];
            }
            sum = sum / maxTime;
            std::cout<<" mean = "<<sum;
            mpf_class squareExpectation = 0;
            for(int time = 0; time < maxTime; time++){
                squareExpectation += (data[time][i][0] * data[time][i][0]) / maxTime;
            }
            squareExpectation -= (sum*sum);
            std::cout<<", var = "<<squareExpectation<<std::endl;
        }
    }

    std::ofstream outputFile;
    outputFile.open("/home/users/r/radhakrb/cardy/data/out.py");
    outputFile<<"data = [";
    for(int time = 0; time < maxTime; time++){
        outputFile<<"[";
        for(int id = 0; id < numWalkers; id++){
            outputFile<<"[";
            for(int d = 0; d < dim; d++){
                outputFile<<data[time][id][d]<<",";
            }
            outputFile.seekp(-1, std::ios_base::cur);
            outputFile<<"],";
        }
        outputFile.seekp(-1, std::ios_base::cur);
        outputFile<<"],";
    }
    outputFile.seekp(-1, std::ios_base::cur);
    outputFile<<"]";



    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    return 0;
}