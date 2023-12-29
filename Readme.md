# TP Micro-aevol2

## Members
Julien PREIN 
Jorge KORGUT Junior

## Project architecture

The project contains 2 majors optimizations and they can be found inside their respectives branches :

* **final-unchanged** : This branch, contains the first version without any modification. This was the starting point of our project.

* **final-custombitset** : This was the first optimization implementation. The idea was to encode the DNA in a bitset, to gain performance by fastly comparing the wanted sequences.

* **final-cuda** : After gaining performance with bitset, the main goal of this step was to convert all the particularities from bitset to the GPU. The propretary language CUDA was used.

## Launching the project

For being able to launch our project, you will need CUDA nvcc compiler installed and the boost library installed inside boost folder.
**Cuda** : ```https://docs.nvidia.com/cuda/cuda-installation-guide-linux```  

**Boost** : ```https://www.boost.org/doc/libs/1_84_0/more/getting_started```  
Extract the files from boost following the example.  
Ex : ```boost/boost_1_83_0```  




