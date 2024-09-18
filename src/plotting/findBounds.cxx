#include<iostream>
#include<fstream>
#include<algorithm>
#include<string>
#include<cmath>

int main(int argc, char **argv){
    float leftBound = std::numeric_limits<float>::infinity();
    float rightBound = -std::numeric_limits<float>::infinity();

    std::string line;
    std::ifstream file(argv[1]);
    if (file.is_open()) {
        while (getline(file, line)) {
            leftBound = std::min(stof(line), leftBound);
            rightBound = std::max(stof(line), rightBound);
        }
        file.close();
    }

    std::cout<<"Minimum is "<<leftBound<<", maximum is "<<rightBound<<std::endl;

    return 0;
}