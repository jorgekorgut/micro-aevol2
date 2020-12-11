add_library(cuda_micro_aevol
        cuExpManager.cu
        cuIndividual.cu
        cuProtein.cu
        misc_functions.cu
        RandService.cu)

set_target_properties(cuda_micro_aevol PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_include_directories(cuda_micro_aevol PUBLIC ../)
target_compile_definitions(cuda_micro_aevol PUBLIC USE_CUDA)