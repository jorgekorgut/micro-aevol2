add_library(cuda_micro_aevol
        Algorithms.cu
        ExpManager_CUDA.cu
        Individual_CUDA.cu
        Protein_CUDA.cu
        misc_functions.cu)

target_include_directories(cuda_micro_aevol PUBLIC ../)