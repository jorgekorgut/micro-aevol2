cmake_minimum_required(VERSION 3.16.3)
project(pdc_mini_aevol)

set(CMAKE_CXX_STANDARD 14)

if ( DO_TRACES )
    add_definitions(-DTRACES)
    message( STATUS "Traces are activated" )
endif ( DO_TRACES )

if ( USE_CUDA )
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        set(CMAKE_CUDA_STANDARD 14)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        set(CMAKE_CUDA_FLAGS "-arch=compute_60")
    else()
        message(FATAL_ERROR "No CUDA support found")
    endif()
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

if ( USE_CUDA )
    add_definitions(-DUSE_CUDA)
    add_subdirectory(cuda)
    add_executable(micro_aevol_gpu main.cpp)
    # nvToolsExt for enhanced profiling (ad-hoc chunks)
    target_link_libraries(micro_aevol_gpu PUBLIC cuda_micro_aevol micro_aevol nvToolsExt)
else ()
    add_executable(micro_aevol_cpu main.cpp)
    target_link_libraries(micro_aevol_cpu micro_aevol)
endif ()