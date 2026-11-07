cmake_minimum_required(VERSION 4.1)

project(numfort_client LANGUAGES Fortran)

find_package(numfort CONFIG REQUIRED)

add_executable(client client.f90)
target_link_libraries(client PRIVATE numfort::numfort)
