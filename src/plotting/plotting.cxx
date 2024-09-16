#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<cmath>
#include<cassert>

int main(int argc, char **argv){
    float leftBound = atof(argv[2]);
    float rightBound = atof(argv[3]);
    int bins = atof(argv[4]);

    std::ofstream tfile;
    tfile.open(argv[6]);

    float resolution = (rightBound - leftBound)/static_cast<float>(bins);

    std::vector<int> counts(bins, 0);
    // for(auto x: counts){
    //     tfile<<x<<std::endl;
    // }

    float data;
    int loc;

    std::string line;
    std::ifstream file(argv[1]);

    if (file.is_open()) {
        while (getline(file, line)) {
            data = stof(line);
            loc = static_cast<int>(floor((data - leftBound)/resolution));
            // tfile<<loc<<", "<<counts[loc]<<std::endl;
            assert(loc>=0 && loc<=bins);
            counts[loc] += 1;
        }
        file.close();
    }

    std::ofstream outputFile;
    outputFile.open(argv[5]);

    for(int i = 0; i < bins; i++){
        outputFile<<i<<"\t"<<counts[i]<<"\n";
    }

    return 0;
}