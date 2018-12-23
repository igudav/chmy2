cmake_minimum_required(VERSION 3.10)
project(chmy2 C)

set(CMAKE_C_STANDARD 11)

add_executable(chmy2 main.c)

target_link_libraries(chmy2 m)