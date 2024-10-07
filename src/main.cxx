#include "main.hpp"

#define PREC 500

simulationParams thisSimulation = {
    .dim = 100,
    .maxTime = 700000,
    .numWalkers = 10,
    
    .low = 200,
    .high = 300,
    .step_size = 0.1,
    .boundary = 500,

    .ccharge = 10,
    .beta = 78.9568,
};

int main(int argc, char **argv){
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    mpf_set_default_prec(PREC);
    mpfr_set_default_prec(PREC);

    std::string dir = argv[1];
    std::filesystem::create_directory(dir);

    std::ofstream parametersFile;
    parametersFile.open(dir + "/parameters.txt");
    thisSimulation.printParameters(parametersFile);
    parametersFile.close();

    std::string parquetFile = dir + "/data.parquet";

    std::vector<std::vector<float>> batchData(
    thisSimulation.maxTime,
    std::vector<float>(thisSimulation.dim*thisSimulation.numWalkers, 0));

    int ratioCount = 0;

    #pragma omp parallel
    {
        #pragma omp for
        for(int id = 0; id < thisSimulation.numWalkers; id++){
            randomWalker walker(thisSimulation, id, read_urandom());
            walker.writeData(batchData, 0);
            for(int t = 1; t < thisSimulation.maxTime; t++){
                walker.cardyEvolve(thisSimulation);
                walker.writeData(batchData, t);
            }

        }
    }
    
    writeBatchToParquet(batchData, parquetFile);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    return 0;
}