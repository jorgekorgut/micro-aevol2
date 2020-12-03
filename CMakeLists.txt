cmake_minimum_required(VERSION 3.16.3)
project(pdc_mini_aevol)

set(CMAKE_CXX_STANDARD 14)

if ( DO_TRACES )
    add_definitions(-DTRACES)
    message( STATUS "Traces are activated" )
endif ( DO_TRACES )

if ( USE_CUDA )
    find_package(CUDA REQUIRED)
    enable_language(CUDA)
    set(CUDA_NVCC_FLAGS "-std=c++11 -arch=compute_60 -lineinfo")
    cuda_compile(cuda_aevol Algorithms.cu)
endif ( USE_CUDA )

if ( USE_OMP )
    find_package(OpenMP REQUIRED)
endif ( USE_OMP )

find_package(ZLIB REQUIRED)

add_library(micro_aevol
        Abstract_ExpManager.cpp
        ExpManager.cpp
        AeTime.cpp
        DnaMutator.cpp
        MutationEvent.cpp
        Organism.cpp
        Stats.cpp
        Threefry.cpp
        Dna.cpp)

target_link_libraries(micro_aevol PUBLIC ZLIB::ZLIB)

if ( OPENMP_FOUND )
    add_definitions(-DUSE_OMP)
    target_link_libraries(micro_aevol PUBLIC OpenMP::OpenMP_CXX)
endif ()


if ( CUDA_FOUND )
    add_definitions(-DUSE_CUDA)
    cuda_add_executable(micro_aevol_gpu main.cpp ${cuda_aevol})
    # nvToolsExt for enhanced profiling (ad-hoc chunks)
    target_link_libraries(micro_aevol_gpu micro_aevol nvToolsExt)
else ( CUDA_FOUND )
    add_executable(micro_aevol_cpu main.cpp)
    target_link_libraries(micro_aevol_cpu micro_aevol)
endif ( CUDA_FOUND )
