#ifndef UTILS_H
#define UTILS_H

#include<fstream>

#include<arrow/api.h>
#include<arrow/io/api.h>
#include<parquet/arrow/writer.h>

void writeBatchToParquet(const std::vector<std::vector<float>>& batchData, const std::string& filename);

unsigned long long int read_urandom();

#endif // UTILS_H
