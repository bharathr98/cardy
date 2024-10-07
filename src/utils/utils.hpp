#ifndef UTILS_H
#define UTILS_H

#pragma once
#include<fstream>
#include<iostream>
#include<string>

#include<gmpxx.h>
#include<mpfr.h>
#include<arrow/api.h>
#include<arrow/io/api.h>
#include<parquet/arrow/writer.h>

void log(std::string message);
struct simulationParams{
    int dim, maxTime, numWalkers;

    mpf_class low, high, step_size, boundary;

    float ccharge, beta;

    void printParameters(std::ofstream& parametersFile){
        parametersFile<<"Precision = "<<mpfr_get_default_prec()<<"\n"
                  <<"Intialisation low bound = "<<low<<"\n"
                  <<"Intialisation high bound = "<<high<<"\n"
                  <<"Step size = "<<step_size<<"\n"
                  <<"Central charge = "<<ccharge<<"\n"
                  <<"Temperature = "<<beta<<"\n"
                  <<"Time steps = "<<maxTime<<"\n"
                  <<"Number of Walkers = "<<numWalkers<<"\n"
                  <<"Dimensions = "<<dim<<"\n";
    }

};

void writeBatchToParquet(const std::vector<std::vector<float>>& batchData, const std::string& filename);

uint64_t read_urandom();

#endif // UTILS_H
