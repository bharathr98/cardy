#ifndef UTILS_H
#define UTILS_H

#pragma once
#include<fstream>

#include<gmpxx.h>
#include<arrow/api.h>
#include<arrow/io/api.h>
#include<parquet/arrow/writer.h>

struct simulationParams{
    int dim, maxTime;

    mpf_class low, high, step_size, boundary;

    float ccharge, beta;

};

void writeBatchToParquet(const std::vector<std::vector<float>>& batchData, const std::string& filename);

uint64_t read_urandom();

#endif // UTILS_H
