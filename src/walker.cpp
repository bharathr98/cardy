#include"boost/multi_array.hpp"
#include<gmp.h>
#include<gmpxx.h>
#include"utils/utils.hpp"
#include"potentials/potentials.hpp"

typedef boost::multi_array<double, 3> array_type;
  array_type A(boost::extents[3][4][2]);

typedef boost::multi_array<mpf_class, 1> step_type; // 1 here is the depth.

class randomWalker{
    private: 
        int maxTime;
        int dim;
        step_type currentStep;
        step_type nextStep;
        mpf_class step_size;
        gmp_randclass rand;

    public:
        randomWalker(simulationParams& _simulationParams, uint64_t seed)
            : currentStep(boost::extents[dim]),
            nextStep(boost::extents[dim]),
            rand(gmp_randinit_default) 
        {
            // Set the parameters
            dim = _simulationParams.dim;
            maxTime = _simulationParams.maxTime;
            step_size = _simulationParams.step_size;

            // Initialise and seed the random number generator
            rand.seed(seed);

            // Populate the currentStep with a randomly chosen point 
            for(int i = 0; i < dim; i++){
                currentStep[i] = (_simulationParams.high - _simulationParams.low)*rand.get_f() + _simulationParams.low;
            }
        }

        void cardyEvolve(simulationParams& _simulationParams){
            for(int d = 0; d < dim; d++){
                // nextStep[d] = abs(prev[i][d] + (2*step_size)*rand.get_f() - step_size);
                nextStep[d] = abs(currentStep[d] + (2*_simulationParams.step_size)*rand.get_f() - _simulationParams.step_size);
                if (nextStep[d] > _simulationParams.boundary){
                    nextStep[d] = 2*_simulationParams.boundary - nextStep[d];
                }
            }
            mpf_class ratioPotential = cardy(_simulationParams.ccharge, _simulationParams.beta, nextStep)/cardy(_simulationParams.ccharge, _simulationParams.beta, currentStep);
            mpf_class decisionToss = rand.get_f();

            if (decisionToss < ratioPotential){
                    currentStep = nextStep;
            }
        }

};