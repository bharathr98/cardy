#include<iostream>
#include<vector>
#include<cmath>

float vanderMonde(std::vector<float> evals){
    int len = evals.size();
    float vm = 1;
    for (int i = 0; i < len; i++){
        for (int j = i + 1; j< len; j++){
            vm *= std::pow((evals[i] - evals[j]),2);
        }
    }
    return vm;
}

int main(){
    std::vector<float> evals = {1,1,2};
    std::cout<<vanderMonde(evals);
    return 0;
}