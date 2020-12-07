add_library(cuda_micro_aevol
        cuExpManager.cu
        cuIndividual.cu
        cuProtein.cu
        misc_functions.cu)

target_include_directories(cuda_micro_aevol PUBLIC ../)