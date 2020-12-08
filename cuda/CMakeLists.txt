add_library(cuda_micro_aevol
        cuExpManager.cu
        cuIndividual.cu
        cuProtein.cu
        misc_functions.cu)

set_target_properties(cuda_micro_aevol PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_include_directories(cuda_micro_aevol PUBLIC ../)