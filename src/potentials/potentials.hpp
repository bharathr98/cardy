#pragma once

#include"boost/multi_array.hpp"

#include<gmp.h>
#include<gmpxx.h>
#include<mpfr.h>

struct cardyParams;

mpf_class vanderMonde(boost::multi_array<mpf_class, 1>& evals);

mpf_class cardyError(float ccharge, float beta, boost::multi_array<mpf_class, 1>& evals);

mpf_class expV(float ccharge, float beta, int dim, boost::multi_array<mpf_class, 1>& evals);

mpf_class expV_int_trap(float ccharge, float betamin, float betamax, unsigned int dim, boost::multi_array<mpf_class, 1>& evals, int num_steps);

mpf_class gaussian(float betaReg, boost::multi_array<mpf_class, 1> evals);

mpf_class gaussianwall(float betaReg, const boost::multi_array<mpf_class, 1>& evals, float Ewall);

mpf_class wigner(float betaReg, boost::multi_array<mpf_class, 1> evals);

mpf_class cardy(float ccharge, float beta, int dim, boost::multi_array<mpf_class, 1>& evals);

mpf_class dampedCardy(float ccharge, float beta, float betaReg, boost::multi_array<mpf_class, 1> evals);

mpf_class CardyGwall(float ccharge, float beta, int dim, boost::multi_array<mpf_class, 1> evals);

mpf_class Cardy_int_trap(float ccharge, float betamin, float betamax, int dim, boost::multi_array<mpf_class, 1> evals);

mpf_class expV_GK(float ccharge, int order, unsigned int dim, boost::multi_array<mpf_class, 1>& evals);

mpf_class Cardy_int_GK(float ccharge, int order, int dim, float Ewall, boost::multi_array<mpf_class, 1> evals);