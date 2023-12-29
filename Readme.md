# TP Micro-aevol2

## Members
Julien PREIN  
Jorge KORGUT Junior

## Project architecture

The project contains 2 majors optimizations and they can be found inside their respectives branches :

* **final-unchanged** : This branch, contains the first version without any modification. This was the starting point of our project.

* **final-custombitset** : This was the first optimization implementation. The idea was to encode the DNA in a bitset, to gain performance by fastly comparing the wanted sequences.

* **final-cuda** : After gaining performance with bitset, the main goal of this step was to convert all the particularities from bitset to the GPU. The propretary language CUDA was used.

### Dependencies

For being able to launch our project, you will need CUDA nvcc compiler installed and the boost library installed inside boost folder.  

**Cuda** : [CUDA installation guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux)

**Boost** : [BOOST Getting started](https://www.boost.org/doc/libs/1_84_0/more/getting_started)  

Extract and install the files from boost following the following folder hierarchic example. The boost folder can be found in the root of the project. 
Ex : ```boost/boost_1_83_0```  

#### :warning: Known errors
* 1
```
error: parameter packs not expanded with ‘...’
```  

>The error comes probably because you have installed cuda-toolkit from the unix package manager. The solution I found to this error is to uninstall every package related to cuda and then reinstall it following the exact instructions from [CUDA installation guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux)  

* 2
```
```
> The reason for this error is the missing CUDA path directory for **nvcc** and the CUDA library. To fix it just add this lines into your ~/.bashrc 
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




