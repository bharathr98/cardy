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


#define PREC 3000


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

mpfr_rnd_t rnd = mpfr_get_default_rounding_mode();

mpf_class expV(float ccharge, float beta, boost::multi_array<mpf_class, 1> evals){
    mpfr_t ccharger, betar, betapr;
    mpfr_init(ccharger);
    mpfr_set_flt(ccharger, ccharge, rnd);
    mpfr_init(betar);
    mpfr_set_flt(betar, beta, rnd);
    mpfr_init(betapr);
    mpfr_set_flt(betapr, 4*M_PI*M_PI/beta, rnd);
    
    // This implements Z = sum(exp(lambda beta)) and Zt = sum(exp(lambda beta,))
    mpfr_t Z;
    mpfr_init_set_str(Z, "0", 10, rnd); // Initialise Z to 0 in base 10
    mpfr_t temp;
    mpfr_init(temp);

    mpfr_t Zt;
    mpfr_init_set_str(Zt, "0", 10, rnd); // Initialise Z to 0 in base 10
    mpfr_t tempt;
    mpfr_init(tempt);

    for (auto x: evals){
        mpfr_set_f(temp, x.get_mpf_t(), rnd);
        mpfr_mul(temp, temp, betar, rnd);
        mpfr_exp(temp, temp, rnd);
        mpfr_add(Z, Z, temp, rnd);

        mpfr_set_f(tempt, x.get_mpf_t(), rnd);
        mpfr_mul(tempt, tempt, betapr, rnd);
        mpfr_exp(tempt, tempt, rnd);
        mpfr_add(Zt, Zt, tempt, rnd);
    }
    mpfr_printf ("Z is %.60Rf\n", Z);
    mpfr_printf ("Zt is %.60Rf\n", Zt);
    

    // This implements exp(Z - Zt)
    mpfr_t potential;
    mpfr_init(potential);
    mpfr_set(potential, Z, rnd);
    mpfr_sub(potential, potential, Zt, rnd);
    mpfr_exp(potential, potential, rnd);


    mpf_t returnValue;
    mpf_init(returnValue);
    mpfr_get_f(returnValue, potential, rnd);
    return mpf_class(returnValue);
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
    mpf_set_default_prec(PREC);
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    gmp_randclass rand (gmp_randinit_default);

    mpf_class low, high;
    low = 0;
    high = 10;

    // Create a 3D array that is time x id x dim
    int time = 3;
    int id = 4;
    int dim = 100;
    typedef boost::multi_array<mpf_class, 3> array_type; // 3 here is the depth.
    boost::array<array_type::index, 3> shape = {{ time, id, dim }}; // 3 here seems to do nothing as long as it is greater than depth. 
    array_type data(shape);
    
    // Randomly initialise the first time step
    // for(int i = 0; i < id; i++){
    //     for(int j = 0; j < dim; j++){
    //         rand.seed(read_urandom());
    //         data[0][i][j] = (high - low)*rand.get_f() + low;
    //     }
    // }
    for(int i = 0; i < id; i++){
        for(int j = 0; j < dim; j++){
            rand.seed(read_urandom());
            data[0][i][j] = (high - low)*rand.get_f() + low;
        }
    }
    std::cout<<data[0][0][0]<<std::endl;
    std::cout<<expV(1, 6.3, data[0][0]);


    // for(int i = 0; i < dimensions; i++){
    // rand.seed(read_urandom());
    // vec[i] = (high - low)*rand.get_f() + low;
    // }
    // std::cout<<vec[5]<<"\n";
    // std::cout<<vanderMonde(vec)<<"\n";
    // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    // std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    return 0;
}