# Mini-Aevol : A mini-application based on the Aevol simulator

A reduced version (from a model and implementation point of view) of Aevol.

DO NOT USE IT TO SIMULATE BIOLOGICAL RESULTS ! See [http://www.aevol.fr](http://www.aevol.fr) for that !

It must be used only to test HPC optimization of the code (parallel, vector, porting to new architecture...).

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

What things you need to install the software and how to install them.
First, you will need to install zlib (and its headers).

```
sudo apt install zlib1g-dev
```

### Compilation

The compilation is straightforward

```
cmake .
make
```

It will produced the executable pdc_micro_aevol.

If CUDA toolkit is available on your system, you can build the software with GPU support

```
cmake . -DUSE_CUDA=on
make
```

It will produced the executable pdc_micro_aevol_gpu.

## Running a simulation

A help (-e or --help) is given to explain the different parameters.

Basically, you must create a directory to store the simulation files (backup/checkpointing and stats files) and then run the simulation
```
mkdir simulation_example_1
cd simulation_example_1
../pdc_micro_aevol
```

You can also resume a simulation from a backup/checkpointing file (for example, resuming from the generation 1000):
```
cd simulation_example_1
../pdc_mini_aevol -r 1000
```

## Model and Implementation

More details about the model and its implementation are given at : 

## Authors

* **Jonathan Rouzaud-Cornabas** - *Initial work*
* **Laurent Turpin**

For the authors of Aevol software, see [http://www.aevol.fr](http://www.aevol.fr)

## License

This project is licensed under the GPLv2 License
