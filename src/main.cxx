#include "main.hpp"

#define PREC 500

namespace po = boost::program_options;

simulationParams thisSimulation;

int main(int argc, char **argv){

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "Produces this help message")
        ("dim", po::value<int>(&thisSimulation.dim)->default_value(100), "Dimensions of the walker")
        ("iterations", po::value<int>(&thisSimulation.maxTime)->default_value(100'000), "Total number of iterations")
        ("initLow", po::value<float>(&thisSimulation.low_f)->default_value(10), "Initialisation Lower Bound")
        ("initHigh", po::value<float>(&thisSimulation.high_f)->default_value(100), "Initialisation Higher Bound")
        ("stepSize", po::value<float>(&thisSimulation.step_size_f)->default_value(0.02), "Walker step size")
        ("boundary", po::value<float>(&thisSimulation.boundary)->default_value(150), "Value at which the Gaussian wall activates")
        ("ccharge", po::value<float>(&thisSimulation.ccharge)->default_value(2), "Central Charge")
        ("gkorder", po::value<int>(&thisSimulation.gk_order)->default_value(50), "Order of the Gauss Kronrod integral. Current allowed values are {17, 18, 50}")
        ("dir", po::value<std::string>(), "Output directory")
    ;

    po::positional_options_description p;
    p.add("dir", -1);

    po::variables_map params;
    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), params);
    po::notify(params);

    // po::variables_map params;
    // po::store(po::parse_command_line(argc, argv, desc), params);
    // po::notify(params);    

    if (params.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    thisSimulation.set_mpf();
    thisSimulation.numWalkers = omp_get_max_threads() - 1;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    mpf_set_default_prec(PREC);
    mpfr_set_default_prec(PREC);

    std::string dir = params["dir"].as<std::string>();
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

    // #pragma omp parallel
    // {   
    //     #pragma omp for schedule(dynamic)
    //     for(int id = 0; id < thisSimulation.numWalkers; id++){
    //         randomWalker walker(thisSimulation, id, read_urandom());
    //         walker.writeData(batchData, 0);
    //         for(int t = 1; t < thisSimulation.maxTime; t++){
    //             walker.cardyEvolve(thisSimulation);
    //             walker.writeData(batchData, t);
    //             if (t % 100 == 0){
    //                 std::cout<<t<<std::endl;
    //             }
    //         }

    //     }
    // }
    
    #pragma omp parallel
{   
    int thread_id = omp_get_thread_num();
    
    if (thread_id != 0) {  // Non-master threads (1-5)
        // Each non-master thread takes one walker
        int id = thread_id - 1;  // Walker 0 goes to thread 1, walker 1 to thread 2, etc.
        
        if (id < thisSimulation.numWalkers) {  // Only process if there's a walker to handle
            randomWalker walker(thisSimulation, id, read_urandom());
            walker.writeData(batchData, 0);
            for(int t = 1; t < thisSimulation.maxTime; t++){
                walker.cardyEvolve(thisSimulation);
                walker.writeData(batchData, t);
                if (t % 10000 == 0){
                    std::cout << "Thread " << thread_id << " processing walker " << id 
                             << " at time " << t << std::endl;
                }
            }
        }
    }
}

    writeBatchToParquet(batchData, parquetFile);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    return 0;
}