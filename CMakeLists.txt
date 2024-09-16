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

add_executable(cardy)
target_sources(cardy
                PRIVATE
                src/main.cxx)
target_link_libraries(cardy 
                        PRIVATE
                        ${GMP_LIBRARIES}
                        ${GMPXX_LIBRARIES}
                        ${MPFR_LIBRARIES}
                        ${Boost_LIBRARIES})

add_executable(plotting)
target_sources(plotting
                PRIVATE
                src/plotting/plotting.cxx)

add_executable(findBounds)
target_sources(findBounds
                PRIVATE
                src/plotting/findBounds.cxx)

install(TARGETS cardy plotting findBounds RUNTIME)