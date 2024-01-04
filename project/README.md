# Mini-Aevol : A mini-application based on the Aevol simulator

A reduced version (from a model and implementation point of view) of Aevol.

DO NOT USE IT TO SIMULATE BIOLOGICAL RESULTS ! See [http://www.aevol.fr](http://www.aevol.fr) for that !

It must be used only to test HPC optimization of the code (parallel, vector, porting to new architecture...).

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites
You will need a unix base system. Sorry for Windows users, maybe try with [Cygwin](http://www.cygwin.com/)

You will also need to install zlib (and its headers):
+ On Debian, Ubuntu, Mint, ... (apt): `zlib1g-dev`
+ On Arch (pacman): `zlib`
+ On Fedora (dnf, rpm): `zlib-devel`

To use cuda implementation, you will need to install the cuda-toolkit:
+ On all distributions: `cuda`

### Compilation

The compilation is straightforward
```
mkdir build
cd build
cmake ..
make
```
It will produced the executable micro_aevol_cpu.

If CUDA toolkit is available on your system, you can build the software with GPU support
```
cmake .. -DUSE_CUDA=on
make
```
It will produced the executable micro_aevol_gpu.

## Running a simulation

A help is given to explain the different parameters when using option `-H` or `--help`.

Basically, you must create a directory to store the simulation files (backup/checkpointing and stats files) and then run the simulation
```
mkdir simulation_example_1
cd simulation_example_1
PATH/TO/micro_aevol_cpu
```

You can also resume a simulation from a backup/checkpointing file (for example, resuming from the generation 1000):
```
cd simulation_example_1
PATH/TO/micro_aevol_cpu -r 1000
```

## Model and Implementation

These [slides](/presentation/slides.pdf) give a short presentation of the model and the purpose of this project can

## Authors

* **Jonathan Rouzaud-Cornabas** - *Initial work*
* **Laurent Turpin**

For the authors of Aevol software, see [http://www.aevol.fr](http://www.aevol.fr)

## License

This project is licensed under the GPLv2 License
