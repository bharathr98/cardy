cmake_minimum_required(VERSION 3.16)

project(cardy LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")
message(${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_BUILD_RPATH "${CMAKE_CURRENT_SOURCE_DIR}/build")

add_compile_options(-lgmp)

find_package(GMP REQUIRED)
if(GMP_FOUND AND GMP_INCLUDE_DIR)
    include_directories(${GMP_INCLUDE_DIR})
endif()

find_package(GMPXX REQUIRED)
if(GMPXX_FOUND AND GMPXX_INCLUDE_DIR)
    include_directories(${GMPXX_INCLUDE_DIR})
endif()

find_package(MPFR REQUIRED)
if(MPFR_FOUND AND MPFR_INCLUDE_DIR)
    include_directories(${MPFR_INCLUDE_DIR})
endif()

find_package(Boost REQUIRED)
message(STATUS "Boost version: ${Boost_VERSION}")
if(Boost_FOUND AND Boost_INCLUDE_DIR)
    include_directories(${Boost_INCLUDE_DIR})
endif()

find_package(Arrow REQUIRED)
if(Arrow_FOUND AND Arrow_INCLUDE_DIR)
    include_directories(${Arrow_INCLUDE_DIR})
endif()

find_package(Parquet REQUIRED)
if(Parquet_FOUND AND Parquet_INCLUDE_DIR)
    include_directories(${Parquet_INCLUDE_DIR})
endif()

add_executable(cardy)
target_sources(cardy
                PRIVATE
                src/main.cxx
                src/utils/read_urandom.cpp
                src/utils/parquet_writer.cpp
                src/potentials/potentials.cpp
                src/walker.cpp)
target_link_libraries(cardy 
                        PRIVATE
                        ${GMP_LIBRARIES}
                        ${GMPXX_LIBRARIES}
                        ${MPFR_LIBRARIES}
                        ${Boost_LIBRARIES}
                        ${Arrow_LIBRARIES}
                        ${Parquet_LIBRARIES}
                        Arrow::arrow_shared
                        Parquet::parquet_shared)

add_executable(plotting)
target_sources(plotting
                PRIVATE
                src/plotting/plotting.cxx)

add_executable(findBounds)
target_sources(findBounds
                PRIVATE
                src/plotting/findBounds.cxx)

add_executable(bench)
target_sources(bench
                PRIVATE
                bench/urandom.cpp
                src/utils/read_urandom.cpp)
target_link_libraries(bench 
                        PRIVATE
                        ${GMP_LIBRARIES}
                        ${GMPXX_LIBRARIES}
                        ${Arrow_LIBRARIES}
                        ${Parquet_LIBRARIES}
                        Arrow::arrow_shared
                        Parquet::parquet_shared)

install(TARGETS cardy plotting findBounds RUNTIME)