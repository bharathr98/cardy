#include"boost/multi_array.hpp"
#include<gmpxx.h>
#include"utils/utils.hpp"
#include"potentials/potentials.hpp"

typedef boost::multi_array<mpf_class, 1> step_type; // 1 here is the depth.

class randomWalker{
    private: 
        int maxTime;
        int dim;
        int id;
        step_type currentStep;
        step_type nextStep;
        mpf_class step_size;
        gmp_randclass rand;
    public:
        randomWalker(simulationParams& _simulationParams, int _id, uint64_t seed);
        void cardyEvolve(simulationParams& _simulationParams);
        void writeData(std::vector<std::vector<float>>& batchData, int time_step);
};