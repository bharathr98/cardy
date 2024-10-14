#include"boost/multi_array.hpp"
#include<gmp.h>
#include<gmpxx.h>
#include<iostream>
#include"utils/utils.hpp"
#include"potentials/potentials.hpp"
#include"walker.hpp"

randomWalker::randomWalker(simulationParams& _simulationParams, int _id, uint64_t seed)
            : currentStep(boost::extents[_simulationParams.dim]),
            nextStep(boost::extents[_simulationParams.dim]),
            rand(gmp_randinit_default) 
        {
            // Set the parameters
            dim = _simulationParams.dim;
            maxTime = _simulationParams.maxTime;
            step_size = _simulationParams.step_size;
            id = _id;

            // Initialise and seed the random number generator
            rand.seed(seed);

            // Populate the currentStep with a randomly chosen point 
            for(int i = 0; i < dim; i++){
                currentStep[i] = (_simulationParams.high - _simulationParams.low)*rand.get_f() + _simulationParams.low;
            }
        }

void randomWalker::cardyEvolve(simulationParams& _simulationParams){
            for(int d = 0; d < dim; d++){
                nextStep[d] = abs(currentStep[d] + (2*_simulationParams.step_size)*rand.get_f() - _simulationParams.step_size);
//                if (nextStep[d] > _simulationParams.boundary){
//                    nextStep[d] = 2*_simulationParams.boundary - nextStep[d];
//                }
            }
            mpf_class denominator = Cardy_int_gk(_simulationParams.ccharge, 6.29, _simulationParams.beta, dim, currentStep);
            if (denominator == 0) {
                std::cerr << "Error: cardy returned zero for currentStep, division by zero!" << std::endl;
                return; // or handle it in some other way
            }
            mpf_class ratioPotential = Cardy_int_gk(_simulationParams.ccharge, 6.29,_simulationParams.beta, dim, nextStep) / denominator;
            // mpf_class ratioPotential = cardy(_simulationParams.ccharge, _simulationParams.beta, nextStep)/cardy(_simulationParams.ccharge, _simulationParams.beta, currentStep);
            mpf_class decisionToss = rand.get_f();

            if (decisionToss < ratioPotential){
                    currentStep = nextStep;
            }
        }
    
void randomWalker::writeData(std::vector<std::vector<float>>& batchData, int time_step){
    for (int d = 0; d < dim; d++){
        batchData[time_step][d + dim*id] = (float)(currentStep[d]).get_d();
    }
}