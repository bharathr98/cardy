#include "utils.hpp"

void writeBatchToParquet(const std::vector<std::vector<float>>& batchData, const std::string& filename) {
    // Create Arrow schema with 1000 columns of double type
    std::vector<std::shared_ptr<arrow::Field>> fields;
    int columnNumber = batchData[0].size();
    for (int i = 0; i < columnNumber; ++i) {
        fields.push_back(arrow::field("col_" + std::to_string(i), arrow::float32()));
    }
    std::shared_ptr<arrow::Schema> schema = arrow::schema(fields);

    // Prepare to write multiple rows to the Parquet file
    arrow::MemoryPool* pool = arrow::default_memory_pool();
    std::vector<std::shared_ptr<arrow::Array>> columns(columnNumber); // 1000 columns

    // Initialize builders for each column
    std::vector<std::unique_ptr<arrow::FloatBuilder>> builders(columnNumber);
    for (int i = 0; i < columnNumber; ++i) {
        builders[i] = std::make_unique<arrow::FloatBuilder>(pool);
    }

    // Append the rows of batch data to each column's builder
    for (const auto& row : batchData) {
        for (int col = 0; col < columnNumber; ++col) {
            builders[col]->Append(row[col]);
        }
    }

    // Build the Arrow Arrays from the builders
    for (int i = 0; i < columnNumber; ++i) {
        builders[i]->Finish(&columns[i]);
    }

    // Create a table from the Arrow arrays
    std::shared_ptr<arrow::Table> table = arrow::Table::Make(schema, columns, batchData.size());

    // Open the file output stream for Parquet
    std::shared_ptr<arrow::io::FileOutputStream> outfile;
    PARQUET_ASSIGN_OR_THROW(outfile, arrow::io::FileOutputStream::Open(filename));

    // Write the Arrow table as Parquet
    PARQUET_THROW_NOT_OK(parquet::arrow::WriteTable(*table, pool, outfile, batchData.size()));

    // Close the output stream
    outfile->Close();
}
