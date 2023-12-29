// ***************************************************************************************************************
//
//          Mini-Aevol is a reduced version of Aevol -- An in silico experimental evolution platform
//
// ***************************************************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <https://gitlab.inria.fr/rouzaudc/mini-aevol>
// Web: https://gitlab.inria.fr/rouzaudc/mini-aevol
// E-mail: See <jonathan.rouzaud-cornabas@inria.fr>
// Original Authors : Jonathan Rouzaud-Cornabas
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ***************************************************************************************************************



#include <iostream>
#include <getopt.h>
#include <cstring>

#ifdef USE_CUDA
#include "cuda/cuExpManager.h"
#endif
#include "Abstract_ExpManager.h"
#include "ExpManager.h"

void print_help(char* prog_path) {
    // Get the program file-name in prog_name (strip prog_path of the path)
    char* prog_name; // No new, it will point to somewhere inside prog_path
    if ((prog_name = strrchr(prog_path, '/'))) prog_name++;
    else prog_name = prog_path;

    printf("******************************************************************************\n");
    printf("*                                                                            *\n");
    printf("*                   micro-aevol - Artificial Evolution                       *\n");
    printf("*                                                                            *\n");
    printf("* Aevol is a simulation platform that allows one to let populations of       *\n");
    printf("* digital organisms evolve in different conditions and study experimentally  *\n");
    printf("* the mechanisms responsible for the structuration of the genome and the     *\n");
    printf("* transcriptome.                                                             *\n");
    printf("*                                                                            *\n");
    printf("* This is a mini-version with a cutdown biological model and a lot of        *\n");
    printf("* features missing.                                                          *\n");

    printf("*        IF YOU WANT TO USE IT TO STUDY EVOLUTION, DO NOT DO IT !            *\n");
    printf("*                                                                            *\n");
    printf("******************************************************************************\n");
    printf("\n");
    printf("Usage : %s -H or --help\n", prog_name);
    printf("   or : %s \n", prog_name);
    printf("   or : %s -r GENERATION_STEP\n", prog_name);
    printf("\nOptions\n");
    printf("  -H, --help\tprint this help, then exit\n");
    printf("  -n, --nbsteps\tNumber of generation to run\n");
    printf("  -w, --width WIDTH_SIZE\tWidth of the population grid is WIDTH_SIZE\n");
    printf("  -h, --height HEIGHT_SIZE\tHeight of the population grid is HEIGHT_SIZE\n");
    printf("  -m, --mutation_rate MUTATION_RATE\tMutation rate is set to MUTATION_RATE\n");
    printf("  -g, --genome_size GENOME_SIZE\tGenome at the initial genome is GENOME_SIZE bps\n");
    printf("  -b, --backup_step BACKUP_STEP\tDo a simulation backup/checkpoint every BACKUP_STEP\n");
    printf("  -r, --resume RESUME_STEP\tResume the simulation from the RESUME_STEP generations\n");
    printf("  -s, --seed SEED\tChange the seed for the pseudo random generator\n");
}

int main(int argc, char* argv[]) {

    int nbstep = -1;
    int width = -1;
    int height = -1;
    double mutation_rate = -1;
    int genome_size = -1;
    int resume = -1;
    int backup_step = -1;
    int seed = -1;

    const char * options_list = "Hn:w:h:m:g:b:r:s:";
    static struct option long_options_list[] = {
            // Print help
            { "help",     no_argument,        NULL, 'H' },
            // Number of generations to be run
            { "nsteps",  required_argument,  NULL, 'n' },
            // Width size of the grid
            { "width", required_argument,  NULL, 'w' },
            // Height size of the grid
            { "height", required_argument,  NULL, 'h' },
            // Mutation rate
            { "mutation_rate", required_argument,  NULL, 'm' },
            // Size of the initial genome
            { "genome_size", required_argument,  NULL, 'g' },
            // Resuming from generation X
            { "resume", required_argument,  NULL, 'r' },
            // Backup step
            { "backup_step", required_argument,  NULL, 'b' },
            // Seed
            { "seed", required_argument,  NULL, 's' },
            { 0, 0, 0, 0 }
    };


    // -------------------------------------------------------------------------
    // 3) Get actual values of the command-line options
    // -------------------------------------------------------------------------
    int option;
    while ((option =
                    getopt_long(argc, argv, options_list, long_options_list, NULL))
           != -1) {
        switch (option) {
            case 'H' : {
                print_help(argv[0]);
                exit(EXIT_SUCCESS);
            }
            case 'w' : {
                width = atoi(optarg);
                break;
            }
            case 'h' : {
                height = atoi(optarg);
                break;
            }
            case 'm' : {
                mutation_rate = atof(optarg);
                break;
            }
            case 'g' : {
                genome_size = atoi(optarg);
                break;
            }
            case 'r' : {
                resume = atoi(optarg);
                break;
            }
            case 'b' : {
                backup_step = atoi(optarg);
                break;
            }
            case 's' : {
                seed = atoi(optarg);
                break;
            }
            case 'n' : {
                nbstep = atoi(optarg);
                break;
            }
            default : {
                // An error message is printed in getopt_long, we just need to exit
                printf("Error unknown parameter\n");
                exit(EXIT_FAILURE);
            }
        }
    }

#ifdef USE_CUDA
    printf("Activate CUDA\n");
#endif

    printf("Start ExpManager\n");

    if (resume >= 0) {
        if ((width != -1) || (height != -1)|| (mutation_rate != -1.0) || (genome_size != -1) ||
            (backup_step != -1) || (seed != -1)) {
            printf("Parameter(s) can not change during the simulation (i.e. when resuming a simulation, parameter(s) can not change)\n");
            exit(EXIT_FAILURE);
        }
        if (nbstep == -1) nbstep = 1000;
    } else {
        if (nbstep == -1) nbstep = 1000;
        if (width == -1) width = 32;
        if (height == -1) height = 32;
        if (mutation_rate == -1) mutation_rate = 0.00001;
        if (genome_size == -1) genome_size = 5000;
        if (backup_step == -1) backup_step = 1000;
        if (seed == -1) seed = 566545665;
    }

    Abstract_ExpManager *exp_manager;
    if (resume == -1) {
        exp_manager = new ExpManager(height, width, seed, mutation_rate, genome_size, backup_step);
    } else {
        printf("Resuming...\n");
        exp_manager = new ExpManager(resume);
    }

#ifdef USE_CUDA
    // Not very clean but the goal is to re-use the initialization on the host to transfer data to device
    auto* tmp = dynamic_cast<ExpManager *>(exp_manager);
    exp_manager = new cuExpManager(tmp);
    delete tmp;
#endif

    exp_manager->run_evolution(nbstep);

    delete exp_manager;

    return 0;
}
