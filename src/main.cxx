#include<iostream>
#include<vector>
#include<memory>
#include<gmp.h>
#include<gmpxx.h>
#include<mpfr.h>
#include<chrono>
#include"boost/multi_array.hpp"
#include<string>
#include<filesystem>
#include"utils/utils.hpp"
#include"potentials/potentials.hpp"

#define PREC 500

int main(int argc, char **argv){
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    mpf_set_default_prec(PREC);
    mpfr_set_default_prec(PREC);
    gmp_randclass rand (gmp_randinit_default);

    std::string dir = argv[2];
    std::filesystem::create_directory(dir);

    mpf_class low, high, step_size, boundary;
    low = 200;
    high = 300;
    step_size = argv[1];
    boundary = 500;

    float ccharge = 10;
    float beta = 78.9568;
    float betaReg = 0.0001;

    // Create a 3D array that is maxTime x numWalkers x dim
    int maxTime = 1'001;
    int numWalkers = 10;
    int dim = 100;
    int batchSize = 1000;

    std::ofstream parametersFile;
    parametersFile.open(dir + "/parameters.txt");

    parametersFile<<"Precision = "<<mpfr_get_default_prec()<<"\n"
                  <<"Intialisation low bound = "<<low<<"\n"
                  <<"Intialisation high bound = "<<high<<"\n"
                  <<"Step size = "<<step_size<<"\n"
                  <<"Central charge = "<<ccharge<<"\n"
                  <<"Temperature = "<<beta<<"\n"
                  <<"Gaussian temperature = "<<betaReg<<"\n"
                  <<"Time steps = "<<maxTime<<"\n"
                  <<"Number of Walkers = "<<numWalkers<<"\n"
                  <<"Dimensions = "<<dim<<"\n";

    std::ofstream outputFile;
    outputFile.open(dir + "/data.txt");
    // outputFile<<"step_size = "<<step_size<<"\n data = [";
    std::string parquetFile = dir + "/data.parquet";

    std::vector<std::vector<float>> batchData(
    maxTime,
    std::vector<float>(dim*numWalkers, 0));

    typedef boost::multi_array<mpf_class, 2> array_type; // 2 here is the depth.
    boost::array<array_type::index, 2> shape = {{ numWalkers, dim }}; // 3 here seems to do nothing as long as it is greater than depth. 
    array_type prev(shape);
    
    // Randomly initialise the first time step
    for(int i = 0; i < numWalkers; i++){
        for(int j = 0; j < dim; j++){
            rand.seed(read_urandom());
            prev[i][j] = (high - low)*rand.get_f() + low;
        }
    }
    int ratioCount = 0;
    // Evolve the walks
    typedef boost::multi_array<mpf_class, 1> step_type; // 1 here is the depth.
    boost::array<step_type::index, 1> stepShape = {{ dim }};
    step_type nextStep(stepShape);
    for(int time = 1; time < maxTime; time++){
        for(int i = 0; i < numWalkers; i++){
            for(int d = 0; d < dim; d++){
                // nextStep[d] = abs(prev[i][d] + (2*step_size)*rand.get_f() - step_size);
                nextStep[d] = abs(prev[i][d] + (2*step_size)*rand.get_f() - step_size);
                if (nextStep[d] > boundary){
                    nextStep[d] = 2*boundary - nextStep[d];
                }
            }
            // mpf_class ratioPotential = cardy(ccharge, beta, nextStep)/cardy(ccharge, beta, prev[i]);
            mpf_class ratioPotential = cardy(ccharge, beta, betaReg, nextStep)/cardy(ccharge, beta, betaReg, prev[i]);
            mpf_class decisionToss = rand.get_f();

            if (decisionToss > ratioPotential){
                    continue;
            }
            else{
                prev[i] = nextStep;
            }

            for(int d = 0; d < dim; d++){
                batchData[time][i*dim + d] = (float)(prev[i][d]).get_d();
            }
        }


        for(int id = 0; id < numWalkers; id++){
            for(int d = 0; d < dim; d++){
                outputFile<<prev[id][d]<<"\n";
            }
        }
        // if (time % 500 == 0){
        //     std::cout<<time<<std::endl;
        // }
    }
    writeBatchToParquet(batchData, parquetFile);
    
    outputFile.close();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    return 0;
}