cmake_minimum_required(VERSION 3.10)

project(MUSIC_DOA_OpenMP VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(OpenMP REQUIRED)

add_executable(music_doa main.cpp)

if(OpenMP_CXX_FOUND)
    target_link_libraries(music_doa PUBLIC OpenMP::OpenMP_CXX)
endif()