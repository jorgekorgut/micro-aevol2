# TP Micro-aevol2

## Members
Julien PREIN  
Jorge KORGUT Junior

## Project architecture

The project contains 2 majors optimizations and they can be found inside their respectives branches :

[ :pushpin: **final-unchanged**  ](https://github.com/jorgekorgut/micro-aevol2/tree/final-unchanged) : This branch, contains the first version without any modification. This was the starting point of our project.

[ :pushpin: **final-custombitset**](https://github.com/jorgekorgut/micro-aevol2/tree/final-custombitset) : This was the first optimization implementation. The idea was to encode the DNA in a bitset, to gain performance by fastly comparing the wanted sequences.

[ :pushpin: **final-cuda**](https://github.com/jorgekorgut/micro-aevol2/tree/final-cuda) : After gaining performance with bitset, the main goal of this step was to convert all the particularities from bitset to the GPU. The propretary language CUDA was used.

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

After **[1]** cloning the repository and **[2]** submodules, **[3]** installing boost inside the /boost folder and **[4]** executing the symlink script, you are ready to change branches and start compiling the different versions of our application.

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

### Bitset
The objective of this project is to optimise the run of the MICRO-AEVOL application. This wide interprétation was purpose chosen to allow us to explore different solutions, technologies and hardwares. There are 2 main axes that we encountered and we found worth it to optimise : run time and memory space.
Taking a glance in the initial code and running it on profiler tools (INTEL VTUNE), we remark that a lot of the runtime was spent traversing the DNA sequence. Therefore this had become our first optimisation target.


![Alt text](resources/hot_path_unchanged_vtune.png)  
_Micro-aevol run profiled with Intel VTune showing the application hot paths_  


Diving into the existing code we discover one of the reasons why the application spends a significant amount of time on the DNA sequence. More precisely, extracting the promoter in a determined position.  

* [int Dna::promoter_at(int pos)](https://github.com/jorgekorgut/micro-aevol2/blob/final-unchanged/project/Dna.cpp#L125C1-L125C1)  

In short, the algorithms used iterates over a specific DNA sequence to find patterns. Therefore, using the classical character vector c++ implementation forces the iteration over the array every time we want to reach our objectives.
The first idea that comes to mind and optimize memory and execution time is to replace this structure by a bitset.

A bitset is a fundamental data structure in computer science that is used to represent a fixed-size sequence of binary digits. For our purpose, the number of digits is defined on run-time, then we cannot use the simple one to build our data in, we need to implement a dynamic version of it. With a little bit of research, we step into a dynamic bitset library called BOOST. This allows us to save up some time and discover the meticulous details of its implementation. BOOST implements it as a dynamic array of blocks that are user defined. These blocks are values and represent one unit of storage supporting the same amount of elements as its bit length. For example : a block of type uint64_t has 64 bits of length and supports 64 elements inside of it.


The next step is then to adapt the original MICRO-AEVOL code, replacing the character vector by our dynamic bitset. That structure was used in the program, to represent a sequence of DNA.

![Alt text](resources/lost_performance_naif_dynamic_bitset.png)  
_Micro-aevol run time execution showing lost of performance with boost dynamic bitset_  

The loss of time performance was expected since we simply replaced the structure but we did not use the full power of the bitset by its totality. Even if we did not win the performance we wished for, we saved up a good space in our memory. Indeed, before we used the bit size of a character times the length of the sequence to represent 1 DNA and now we use the ceiling of the length of the sequence divided by the bit size of a block in the dynamic bitset.
To illustrate, if we use a block of uint_64, we economise 8 * n / (n / 64) = 512 times the memory needed to represent the same DNA in the previous version.

Knowing that we can do more, we focus now on modifying the algorithm's logic to work with our new approach. To be precise, we changed the loops of 5 main functions (promoter_at, terminator_at, shine_dal_start, protein_stop, codon_at). The common thing we remark on all these functions is an inner loop to check if a specific pattern is found. Knowing that these patterns are constant and that they are smaller than our max available block bit size, we could heavily optimise these functions. 
To benefit from the optimisations that BOOST's engineers had already think of, we added a new function to their library ([get_subset(size_type pos, size_type len)](https://github.com/druckdev/dynamic_bitset/blob/42271ad9e47df986ec585f8aaa691c34a0819c67/include/boost/dynamic_bitset/dynamic_bitset.hpp#L1401C1) ) that returns a block with the starting bit defined by pos and the last one defined by pos + len - 1. The rest of its bits are 0. This function uses binary operators to quickly retrieve and construct our target block.
With this new helper function, we could adapt our previous functions and surprisingly we had a really nice interface to work with in our use case.
This optimization led our program to execute faster than the previous version and the starting one. 

![Alt text](resources/time_compare_unchanged_dynamic_bitset.png)  

_Run time comparison of the starting version of Micro-Aevol and the one with custom dynamic bitset_  

![Alt text](resources/dynamic_bitset_application_hot_path.png)  

_Micro-aevol custom dynamic bitset version run profiled with Intel VTune showing the application hot paths_  

It is interesting to notice the loss of the relative time spent on the DNA sequence. 
### GPU
