cmake_minimum_required(VERSION 3.23)
project(parallel_multiscale_approximation C)

set(CMAKE_C_STANDARD 99)

add_executable(parallel_multiscale_approximation main.c cuda_test.cu)

target_link_libraries(parallel_multiscale_approximation  m)
target_link_libraries(parallel_multiscale_approximation openblas)