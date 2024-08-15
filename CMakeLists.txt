cmake_minimum_required(VERSION 3.16)

project(cardy LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")
message(${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_BUILD_RPATH "${CMAKE_CURRENT_SOURCE_DIR}/build")

find_package(GMP REQUIRED)
find_package(GMPXX REQUIRED)

add_executable(cardy)
target_sources(cardy
                PRIVATE
                src/main.cxx)
target_link_libraries(cardy 
                        PRIVATE
                        ${GMP_LIBRARIES}
                        ${GMPXX_LIBRARIES})
install(TARGETS cardy RUNTIME)