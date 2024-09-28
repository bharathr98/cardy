#include "../src/main.hpp"
#include "../src/utils/utils.hpp"

int main(){
    gmp_randclass rand (gmp_randinit_default);

    int time = 700'000;
    int dim = 100;
    int numWalkers = 10;

    int iterations = time * dim * numWalkers;

    mpf_class total = 0;
    mpf_class xSquared = 0;
    mpf_class newNumber = 0;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    rand.seed(read_urandom());
    for(int i = 0; i < iterations; i++){
        newNumber = rand.get_f();
        total += newNumber;
        xSquared += (newNumber * newNumber);
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout <<"Time elapsed = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    std::cout<<"Expectation value = "<<total/iterations<<std::endl;
    std::cout<<"Variance value = "<<xSquared/iterations - (total/iterations)*(total/iterations)<<std::endl;


    return 0;
}