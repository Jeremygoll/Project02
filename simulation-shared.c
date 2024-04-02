/**
 * Runs the organism-evolution simulator.
 * 
 * This version runs in parallel using shared memory. Compile with:
 *     gcc -Wall -O3 -march=native -fopenmp simulation-shared.c organism.c -o simulation-shared
 * 
 * ***All added code is marked with a NOTE comment***
 */
#define _GNU_SOURCE

#include <omp.h> // NOTE: added in shared version

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include <unistd.h>

#include "organism.h"


int main(int argc, char* const argv[]) {
    // default command-line options
    size_t iterations = 4*1024*1024;
    unsigned int random_seed = time(NULL);
    double mutation_prob = 0.02;
    bool organism_debug = false;
    size_t num_threads = omp_get_max_threads(); // NOTE: added in shared version

    // parse command line arguments
    int opt;
    while ((opt = getopt(argc, argv, "n:r:m:t:d")) != -1) {
        char* end;
        switch (opt) {
        case 'n': iterations = strtoumax(optarg, &end, 10); break;
        case 'r': random_seed = strtoul(optarg, &end, 10); break;
        case 'm': mutation_prob = atof(optarg); break;
        case 'd': organism_debug = true; break;
        case 't': num_threads = strtoumax(optarg, &end, 10); break; // NOTE: added in shared version
        default:
            fprintf(stderr, "usage: %s [-n num-iterations] [-r random-seed] [-m mutation-prob] [-d] [-t num_threads] input output\n", argv[0]);
            return 1;
        }
    }
    if (optind + 2 != argc || iterations == 0 || num_threads == 0) {
        fprintf(stderr, "usage: %s [-n num-iterations] [-r random-seed] [-m mutation-prob] [-d] [-t num_threads] input output\n", argv[0]);
        return 1;
    }
    const char* input_path = argv[optind];
    const char* output_path = argv[optind+1];

    // load the environment
    environment env;
    if (!environment_load_from_path(&env, input_path, random_seed, num_threads)) { return 1; }
    env.mutation_prob = mutation_prob;
    size_t world_size = env.world_size;
    
    if (organism_debug) {
        organism_print(&env.grid[(world_size*world_size+world_size)/2]); // prints middle grid cell
    }

    // start the timer
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // run the simulation
    for (size_t t = 0; t < iterations; t++) {
        // update all organisms
        // NOTE: run the inner loop in parallel
        #pragma omp parallel for default(none) num_threads(num_threads) firstprivate(world_size) shared(env) schedule(dynamic, 16)
        for (size_t i = 0; i < world_size*world_size; i++) {
            organism_update(&env.grid[i]);
        }

        // spawn new organisms
        environment_process_spawned_organisms(&env);
        for (size_t i = 0; i < env.spawned[0].size; i++) {
            organism_spawn(&env, &env.spawned[0].data[i]);
        }
        env.spawned[0].size = 0;

        if (organism_debug && t % (1024*1024) == 0) {
            // debugging print; cannot be done in parallel
            environment_print(&env);
        }
    }
    environment_combine_tasks(&env); // NOTE: added in shared version to combine tasks from all threads

    if (organism_debug && iterations % (1024*1024) != 0) { environment_print(&env); }

    // get the elapsed time
    clock_gettime(CLOCK_MONOTONIC, &end);
    double time = end.tv_sec-start.tv_sec+(end.tv_nsec-start.tv_nsec)/1000000000.0;
    printf("Time: %g secs\n", time);

    // save the results
    environment_save_to_path(&env, output_path);

    // print the tasks completed
    size_t* counts = env.tasks_completed[0].counts;
    printf("NOT:     %zu\n", counts[TASK_NOT]);
    printf("NAND:    %zu\n", counts[TASK_NAND]);
    printf("AND:     %zu\n", counts[TASK_AND]);
    printf("OR-NOT:  %zu\n", counts[TASK_OR_NOT]);
    printf("OR:      %zu\n", counts[TASK_OR]);
    printf("AND-NOT: %zu\n", counts[TASK_AND_NOT]);
    printf("NOR:     %zu\n", counts[TASK_NOR]);
    printf("XOR:     %zu\n", counts[TASK_XOR]);
    printf("EQUAL:   %zu\n", counts[TASK_EQUAL]);

    // cleanup
    environment_free(&env);
}
