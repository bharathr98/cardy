#pragma once

#include"boost/multi_array.hpp"

#include<gmp.h>
#include<gmpxx.h>
#include<mpfr.h>

struct cardyParams;

mpf_class vanderMonde(boost::multi_array<mpf_class, 1> evals);

mpf_class expV(float ccharge, float beta, boost::multi_array<mpf_class, 1> evals);

mpf_class gaussian(float betaReg, boost::multi_array<mpf_class, 1> evals);

mpf_class wigner(float betaReg, boost::multi_array<mpf_class, 1> evals);

mpf_class cardy(float ccharge, float beta, boost::multi_array<mpf_class, 1> evals);

mpf_class dampedCardy(float ccharge, float beta, float betaReg, boost::multi_array<mpf_class, 1> evals);