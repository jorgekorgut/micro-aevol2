# TP Micro-aevol2

## Members
Julian PREIN  
Jorge KORGUT Junior

## Project architecture

The project contains 2 majors optimizations that can be found in the following
branches:

[ :pushpin: **final-unchanged**  ](../../tree/final-unchanged) : This branch, contains the first version without any modification. This was the starting point of our project.

[ :pushpin: **final-custombitset**](../../tree/final-custombitset) : This was the first optimization implementation. The idea was to encode the DNA in a bitset, to gain performance by fastly comparing the wanted sequences.

[ :pushpin: **master**](../../tree/master) : After gaining performance with the bitset, the main goal of this step was to convert all the particularities from bitset to the GPU. The propretary language CUDA was used.

### Dependencies

For being able to launch our project, you will need CUDA nvcc compiler installed and the boost library installed inside boost folder.  

**Cuda** : [CUDA installation guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux)

**Boost** : [BOOST Getting started](https://www.boost.org/doc/libs/1_84_0/more/getting_started)  

Extract and install the files from boost following the following folder hierarchic example. The boost folder can be found in the root of the project. 
Ex : ```boost/boost_1_83_0```  

#### :warning: Known errors

* ```error: parameter packs not expanded with ‘...’```  

The error comes probably because you have installed cuda-toolkit from the unix package manager. The solution I found to this error is to uninstall every package related to cuda and then reinstall it following the exact instructions from [CUDA installation guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux)  

* ```Could NOT find CUDA (missing: CUDA_NVCC_EXECUTABLE CUDA_CUDART_LIBRARY)```  

The reason for this error is the missing CUDA path directory for **nvcc** and the CUDA library. To fix it just add this lines into your ~/.bashrc 
```
export PATH=/usr/local/cuda-12.3/bin${PATH:+:${PATH}}
export LD_LIBRARY_PATH=/usr/local/cuda-12.3/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
```

### Launching the project  

After cloning the repository execute `bootstrap.sh` to install and setup our
boost fork dependency.

As mentioned before, every branch is tagged with a prefix _**final-**_ before the optimization type used. Therefore, for launching the project with a specific optimization, change branches to the wished one.

```
git checkout final-#
```  

Then create a new folder inside projects.

```
mkdir project/build
```  

Launch cmake inside the folder you have created.

```
cd project/build
cmake .. -DUSE_CUDA=on
```  

Build the project
```
make -j
```  

To measure the application execution time, you can execute the following unix application.  
```
time project/build/micro_aevol_#
```

## Report  

### Introduction
The objective of this project is to optimize the run of the MICRO-AEVOL application. This wide scope was purpose chosen to allow us to explore different solutions either with hardware and software optimization. In what concerns this project, 2 main axes will be focused on : algorithm run time and memory space.

### Bitset
By analyzing the initial code and running it on profiler tools (INTEL VTUNE), we remark that a lot of the runtime was spent traversing the DNA sequence. Therefore this had become our first optimization target.


![Alt text](resources/hot_path_unchanged_vtune.png)  
_Micro-aevol run profiled with Intel VTune showing the application hot paths_  


Diving into the existing code we discover one of the reasons why the application spends a significant amount of time on the DNA sequence. In summary multiples algorithms used in our application iterates over the DNA to find patterns. For example, in one of them, we need to extract the promoter in a determined position. For doing so, we cover all the DNA iterating position after position and comparing the starting sequence to a reference and constant pattern, saving its distance for later on.

* [int Dna::promoter_at(int pos)](https://github.com/jorgekorgut/micro-aevol2/blob/final-unchanged/project/Dna.cpp#L125C1-L125C1)  

For that reason, using the classical character vector c++ implementation forces the iteration over the array every time we want to reach our objectives. This approach waste a lot of memory and causes unwanted cache misses to fetch data from further and slower memories.

Therefore, the first optimization idea that pops up is to replace this structure by a bitset. But what is a bitset ?

A bitset is a fundamental data structure in computer science that is used to represent a fixed-size sequence of binary digits. 

However, for our purpose, the number of digits is defined on run-time. Thus we are not able to use the simple implementation to build our data in, we are constraint to implement a dynamic version of this data structure.  

With a little bit of research, we step into a dynamic bitset library called BOOST. This allows us to save up some development time and discover the meticulous details of its implementation. BOOST implements the dynamic bitset as a dynamic array of blocks that are user defined. These blocks are values and represent one unit of storage supporting the same amount of elements as in its bit length. For example : a block of type uint64_t has 64 bits of length and supports 64 elements inside of it.

The next step is to adapt the original MICRO-AEVOL code, replacing the character vector by our dynamic bitset. However, this data structure limit our use case to only 2 different genetics codes. If we want to implement the 4 different ones, the dynamic bitset would still be possible, but an adaptation to 2 instead of 1 bits is required.

![Alt text](resources/lost_performance_naif_dynamic_bitset.png)  
_Micro-aevol run time execution showing lost of performance with boost dynamic bitset. This graph was constructed running 20 times each version and then plotting the result._  

The loss of time performance was expected since we simply replaced the structure but we did not use the full potential of the bitset. Even if we did not win the performance we wished for, we saved up a good amount of space in memory. Indeed, before we used the bit length of a character times the length of the sequence to represent 1 DNA and now we use the ceiling of the length of the sequence divided by the bit length of a block in the dynamic bitset.
To illustrate, if we use a block of uint_64, we economize ```8 * n / (n / 64) = 512```. That is **512 times less memory** to represent the same DNA in the previous version.

Knowing the full potential of bitset, the focus is on modifying the algorithm's logic to work with our new approach. To be precise, we changed the loops of 5 main functions :  

 * promoter_at  
 * terminator_at  
 * shine_dal_start  
 * protein_stop  
 * codon_at  

The common thing we remark on all these functions is an inner loop to check if a specific pattern is found. Knowing that these patterns are constant and that they are smaller than our max available block bit size, we could heavily optimize these functions.  

To benefit from the optimizations that BOOST's engineers had already think of, we added a new function to their library :

* [get_subset(size_type pos, size_type len)](https://github.com/druckdev/dynamic_bitset/blob/42271ad9e47df986ec585f8aaa691c34a0819c67/include/boost/dynamic_bitset/dynamic_bitset.hpp#L1401C1) 

This function returns a block with the starting bit defined by pos and the last one defined by pos + len - 1. The rest of its bits are 0. In short, function uses binary operators to quickly retrieve, construct and return our target sequence.  

With this new helper function, we could adapt our previous functions and surprisingly we had a really nice interface to work with.  

This optimization led our program to execute faster than the previous version and the unchanged one. 

![Alt text](resources/time_compare_unchanged_dynamic_bitset.png)  

_Run time comparison of the starting version of Micro-Aevol and the one with custom dynamic bitset. This graph was constructed running 20 times each version and then plotting the result._  

![Alt text](resources/dynamic_bitset_application_hot_path.png)  

_Micro-aevol custom dynamic bitset version run profiled with Intel VTune showing the application hot paths_  

It is interesting to notice the loss of the relative time spent on the DNA sequence and the total speedup obtained using this optimization method. 

### GPU

After optimization the CPU executable we switched over to the GPU. We started by
a very naive replacement of all the `char`, `*int8_t` and `uint` arrays we could
find to a bitset in the form of `block *` (i.e. `uint64_t`). This obviously
broke the functionality of the program.

From here on we started to also slowly port the functionality of every function
using those arrays. This consisted mostly of adjustments in uses of indices and
lengths in loops, but also of more complex parts like introducing specific
atomic operations since different threads might need to access the same bytes
because of the denser memory layout. In this step we also partially reverted
back to `uint` arrays for the arrays holding indices (e.g. `terminator`,
`prot_start`). At the end of this we had working code again.

When comparing the performance we first benchmarked the state of the program
when we got it. For this we just turned on cuda support and increased the C++
version to 20.

```patch
$ git diff origin/dev-master-cuda 80ed9301ddd0 -- CMakeLists.txt
diff --git a/CMakeLists.txt b/CMakeLists.txt
index d55343c24c74..8efb394019db 100644
--- a/CMakeLists.txt
+++ b/CMakeLists.txt
@@ -1,13 +1,11 @@
 cmake_minimum_required(VERSION 3.13)
 project(pdc_mini_aevol)

-set(CMAKE_CXX_STANDARD 20)
-set(CMAKE_BUILD_TYPE Release)
-set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
+set(CMAKE_CXX_STANDARD 14)

 set(DO_TRACES OFF CACHE BOOL "Whether to enable the Traces library")
 set(USE_OMP OFF CACHE BOOL "Whether to enable OpenMP parallelization")
-set(USE_CUDA ON CACHE BOOL "Whether to enable CUDA parallelization")
+set(USE_CUDA OFF CACHE BOOL "Whether to enable CUDA parallelization")

 if ( DO_TRACES )
     add_definitions(-DTRACES)
```

The benchmarked commit can be found on the
[dev-master-cuda](../../tree/dev-master-cuda) branch. This commit would run in
our tests for 28 seconds and our bitset version for 17.2 seconds.

After this we profiled our code with nvprof and got following results:

```
==10520== NVPROF is profiling process 10520, command: ./micro_aevol_gpu
==10520== Profiling application: ./micro_aevol_gpu
==10520== Profiling result:
            Type  Time(%)      Time     Calls       Avg       Min       Max  Name
 GPU activities:   51.09%  8.72178s      1001  8.7131ms  8.5950ms  13.767ms  sparse_meta(unsigned int, cuIndividual*)
                   13.63%  2.32774s         1  2.32774s  2.32774s  2.32774s  clean_population_metadata(unsigned int, cuIndividual*)
                   12.71%  2.17036s      1001  2.1682ms  2.1056ms  3.4977ms  gather_genes(unsigned int, cuIndividual*)
                    8.07%  1.37777s      1001  1.3764ms  1.3479ms  2.3350ms  search_patterns(unsigned int, cuIndividual*)
                    6.86%  1.17070s      1001  1.1695ms  28.512us  2.0323ms  clean_metadata(unsigned int, cuIndividual*)
                    3.46%  590.94ms      1001  590.35us  571.59us  877.67us  find_gene_per_RNA(unsigned int, cuIndividual*)
                    1.73%  294.78ms      1001  294.49us  290.37us  482.63us  compute_fitness(unsigned int, cuIndividual*, double const *)
                    1.26%  214.40ms      1001  214.18us  208.61us  356.20us  compute_phenotype(unsigned int, cuIndividual*)
                    0.47%  80.508ms         2  40.254ms  31.915ms  48.593ms  check_result(unsigned int, cuIndividual*)
                    0.31%  52.315ms      1001  52.263us  50.464us  82.561us  translation(unsigned int, cuIndividual*)
                    0.18%  30.995ms      1001  30.964us  29.600us  41.632us  transcription(unsigned int, cuIndividual*)
                    0.08%  13.131ms      1000  13.130us  12.576us  19.968us  selection(unsigned int, unsigned int, cuIndividual const *, RandService*, int*)
                    0.06%  11.020ms      1000  11.020us  10.464us  16.160us  reproduction(unsigned int, cuIndividual*, int const *, unsigned long const *)
                    0.06%  9.5282ms      1000  9.5280us  8.2880us  17.088us  do_mutation(unsigned int, cuIndividual*, double, RandService*)
                    0.01%  2.4967ms         1  2.4967ms  2.4967ms  2.4967ms  init_device_population(int, unsigned int, unsigned int, int, cuIndividual*, unsigned long*, unsigned char*, unsigned long*, unsigned int*, unsigned long*, unsigned int*, cuRNA*)
                    0.01%  2.4171ms      1001  2.4140us  2.2400us  12.192us  swap_parent_child_genome(unsigned int, cuIndividual*, unsigned long*)
                    0.01%  1.5842ms      1028  1.5410us  1.2800us  3.2640us  [CUDA memcpy DtoH]
                    0.00%  425.78us      1027     414ns     320ns  2.4000us  [CUDA memcpy HtoD]
      API calls:   85.72%  14.7989s     13013  1.1372ms  5.2830us  48.925ms  cudaDeviceSynchronize
                   13.58%  2.34485s      2055  1.1410ms  3.4150us  2.32776s  cudaMemcpy
                    0.31%  53.503ms         1  53.503ms  53.503ms  53.503ms  cudaDeviceSetLimit
                    0.23%  39.463ms     13014  3.0320us  2.0730us  7.5902ms  cudaLaunchKernel
                    0.11%  19.057ms         1  19.057ms  19.057ms  19.057ms  cudaDeviceReset
                    0.03%  4.9337ms        13  379.51us  1.6610us  3.2099ms  cudaFree
                    0.01%  1.7501ms     13013     134ns      85ns  300.55us  cudaGetLastError
                    0.00%  668.90us        13  51.454us  1.7030us  160.60us  cudaMalloc
                    0.00%  203.28us       114  1.7830us     207ns  79.210us  cuDeviceGetAttribute
                    0.00%  17.983us         1  17.983us  17.983us  17.983us  cuDeviceGetName
                    0.00%  11.497us         1  11.497us  11.497us  11.497us  cuDeviceGetPCIBusId
                    0.00%  2.1560us         3     718ns     329ns  1.4750us  cuDeviceGetCount
                    0.00%  1.1050us         2     552ns     224ns     881ns  cuDeviceGet
                    0.00%     532ns         1     532ns     532ns     532ns  cuDeviceTotalMem
                    0.00%     411ns         1     411ns     411ns     411ns  cuModuleGetLoadingMode
                    0.00%     360ns         1     360ns     360ns     360ns  cuDeviceGetUuid
```

As a next step we parallelized `clean_population_metadata` with a very simple
patch to use a thread per individual instead of just one thread doing
everything:

```
diff --git a/project/cuda/cuExpManager.cu b/project/cuda/cuExpManager.cu
index 110ea097d75f..95f4fcffb0e5 100644
--- a/project/cuda/cuExpManager.cu
+++ b/project/cuda/cuExpManager.cu
@@ -344,7 +344,7 @@ void cuExpManager::transfer_to_host() const {
 }
 
 void cuExpManager::device_data_destructor() {
-    clean_population_metadata<<<1, 1>>>(nb_indivs_, device_individuals_);
+    clean_population_metadata<<<ceil((float)nb_indivs_ / 32.0), 32>>>(nb_indivs_, device_individuals_);
     RandService tmp_rand;
     checkCuda(cudaMemcpy(&tmp_rand, rand_service_, sizeof(RandService), cudaMemcpyDeviceToHost));
     checkCuda(cudaFree(tmp_rand.rng_counters));
@@ -568,10 +568,10 @@ __global__ void compute_fitness(uint size, cuIndividual* individuals, const doub
 // Interface Host | Device
 
 __global__ void clean_population_metadata(uint nb_indivs, cuIndividual* individuals) {
-    if (threadIdx.x + blockIdx.x == 0) {
-        for (int i = 0; i < nb_indivs; ++i) {
-            individuals[i].clean_metadata();
-        }
+    // On thread per individual
+    auto indiv_idx = threadIdx.x + blockIdx.x * blockDim.x;
+    if (indiv_idx < nb_indivs) {
+        individuals[indiv_idx].clean_metadata();
     }
 }
```

With this we gained over 2 seconds with a final runtime of 14.9 seconds. Thus we
achieved a speedup of almost 1.9x.

#### What now?

When looking at the output of nvprof one can clearly see that the majority of
the runtime is spent in `sparse_meta`:

```
==10694== NVPROF is profiling process 10694, command: ./micro_aevol_gpu
==10694== Profiling application: ./micro_aevol_gpu
==10694== Profiling result:
            Type  Time(%)      Time     Calls       Avg       Min       Max  Name
 GPU activities:   59.15%  8.70406s      1001  8.6954ms  8.5703ms  13.858ms  sparse_meta(unsigned int, cuIndividual*)
                   14.72%  2.16590s      1001  2.1637ms  2.0989ms  3.4946ms  gather_genes(unsigned int, cuIndividual*)
                    9.32%  1.37108s      1001  1.3697ms  1.3398ms  2.3108ms  search_patterns(unsigned int, cuIndividual*)
                    7.96%  1.17136s      1001  1.1702ms  28.513us  1.9597ms  clean_metadata(unsigned int, cuIndividual*)
                    4.01%  589.85ms      1001  589.26us  572.04us  886.05us  find_gene_per_RNA(unsigned int, cuIndividual*)
                    2.00%  294.39ms      1001  294.10us  289.73us  484.13us  compute_fitness(unsigned int, cuIndividual*, double const *)
                    1.45%  213.97ms      1001  213.76us  208.13us  357.86us  compute_phenotype(unsigned int, cuIndividual*)
                    0.55%  80.416ms         2  40.208ms  31.854ms  48.562ms  check_result(unsigned int, cuIndividual*)
                    0.35%  52.030ms      1001  51.978us  50.113us  82.817us  translation(unsigned int, cuIndividual*)
                    0.21%  30.956ms      1001  30.925us  29.664us  41.728us  transcription(unsigned int, cuIndividual*)
                    0.09%  12.962ms      1000  12.961us  12.352us  20.191us  selection(unsigned int, unsigned int, cuIndividual const *, RandService*, int*)
                    0.08%  11.050ms      1000  11.049us  10.528us  16.352us  reproduction(unsigned int, cuIndividual*, int const *, unsigned long const *)
                    0.06%  9.3883ms      1000  9.3880us  8.2240us  16.512us  do_mutation(unsigned int, cuIndividual*, double, RandService*)
                    0.02%  2.5412ms      1001  2.5380us  2.3670us  12.608us  swap_parent_child_genome(unsigned int, cuIndividual*, unsigned long*)
                    0.02%  2.5118ms         1  2.5118ms  2.5118ms  2.5118ms  init_device_population(int, unsigned int, unsigned int, int, cuIndividual*, unsigned long*, unsigned char*, unsigned long*, unsigned int*, unsigned long*, unsigned int*, cuRNA*)
                    0.01%  1.6073ms      1028  1.5630us  1.2800us  9.5040us  [CUDA memcpy DtoH]
                    0.01%  983.56us         1  983.56us  983.56us  983.56us  clean_population_metadata(unsigned int, cuIndividual*)
                    0.00%  396.70us      1027     386ns     288ns  2.4320us  [CUDA memcpy HtoD]
      API calls:   99.00%  14.7695s     13013  1.1350ms  5.2760us  48.960ms  cudaDeviceSynchronize
                    0.44%  65.061ms         1  65.061ms  65.061ms  65.061ms  cudaDeviceSetLimit
                    0.27%  39.591ms     13014  3.0420us  2.1670us  6.7657ms  cudaLaunchKernel
                    0.13%  18.682ms         1  18.682ms  18.682ms  18.682ms  cudaDeviceReset
                    0.12%  18.269ms      2055  8.8900us  4.0870us  996.16us  cudaMemcpy
                    0.03%  4.9751ms        13  382.70us  1.9580us  3.2105ms  cudaFree
                    0.01%  1.5044ms     13013     115ns      81ns  298.24us  cudaGetLastError
                    0.00%  674.46us        13  51.881us  1.7180us  160.61us  cudaMalloc
                    0.00%  192.25us       114  1.6860us     169ns  79.173us  cuDeviceGetAttribute
                    0.00%  16.932us         1  16.932us  16.932us  16.932us  cuDeviceGetName
                    0.00%  10.792us         1  10.792us  10.792us  10.792us  cuDeviceGetPCIBusId
                    0.00%  2.0680us         3     689ns     295ns  1.4280us  cuDeviceGetCount
                    0.00%     934ns         2     467ns     196ns     738ns  cuDeviceGet
                    0.00%     464ns         1     464ns     464ns     464ns  cuModuleGetLoadingMode
                    0.00%     444ns         1     444ns     444ns     444ns  cuDeviceTotalMem
                    0.00%     296ns         1     296ns     296ns     296ns  cuDeviceGetUuid
```

This is the case since `sparse_meta` only uses 3 threads instead of all of the
available ones. Our plan was to parallelize sparse_meta by combining most of its
work with `search_patterns`. The problem here was that the order of the
resulting lists `terminator_idxs` and `prot_start_idxs` would be different (i.e.
non-deterministic) which changed the end-result. Additionally we couldn't gain
any runtime through our approach. Our implementation can be found on the
[further](../../tree/further) branch. We believe that with more time one could bring
those two functions together and increase the GPU performance even further.

