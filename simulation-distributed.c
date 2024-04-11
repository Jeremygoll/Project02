/**
 * Runs the organism-evolution simulator.
 * 
 * This version runs in parallel using distributed and shared memory. Compile with:
 *     mpicc -Wall -O3 -march=native -fopenmp simulation-distributed.c organism.c -o simulation-distributed
 *
 * ***All added code is marked with a NOTE comment***
 */
#define _GNU_SOURCE

#include <mpi.h> // NOTE: added in distributed version
#include <omp.h> // NOTE: added in shared version

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <inttypes.h>
#include <unistd.h>

#include "organism.h"


// NOTE: added in distributed version
#include <stdint.h>
#include <limits.h>
#if SIZE_MAX == UCHAR_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "MPI_SIZE_T indeterminate"
#endif

void gather_spawned_organisms_rank_zero(environment* env, int size, MPI_Datatype organism_type) {
    // TODO: every process has its own list of spawned organisms, they need to be merged into a
    // single list in one process before environment_process_spawned_organisms() is called from
    // that process. Then, the list of spawned organisms should be sent back to the other
    // processes so organisms can be spawned in the correct locations. At the end, each process
    // should clear its own list of spawned organisms after the merge by doing
    // env.spawned[i].size = 0; This code should probably be placed in a function.

    // spawn new organisms
    // TODO: need to move all spawned organisms to the primary process before calling environment_process_spawned_organisms()
    environment_combine_spawned_organisms(env); // this sums over all threads in the current process, not across processes
    int num_spawned = env->spawned[0].size;
    int* counts = malloc(sizeof(int) * size); // this is rank 0 only
    MPI_Gather(&num_spawned, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // these are rank 0 only
    int* displacements = (int*)malloc(size * sizeof(int));
    for (int i = 0; i < size; i++) {
        displacements[i] = i == 0 ? 0 : displacements[i - 1] + counts[i - 1];
    }
    int total = displacements[size - 1] + counts[size - 1];

    spawned_org_array* spawned = &env->spawned[0];

    // this is rank 0 only
    spawned_org_array_ensure_capacity(spawned, total);

    MPI_Gatherv(MPI_IN_PLACE, spawned->size, organism_type,
        spawned->data, counts, displacements, organism_type,
        0, MPI_COMM_WORLD);

    environment_process_spawned_organisms(env);
    
    // TODO: Scatter sizes and spawned organisms to all processes
    size_t rows = env->world_size;
    size_t per_process = rows / size;
    size_t remainder = rows % size;
    // count and displacements are being reused from earlier
    for (int i = 0; i < size; i++) {
        counts[i] = per_process + (i < remainder);
        displacements[i] = i * per_process + (i < remainder ? i : remainder);
    }

    MPI_Scatter(counts, 1, MPI_INT, &env->spawned[0].size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(env->spawned[0].data, counts, displacements, organism_type, env->spawned[0].data, env->spawned[0].size, organism_type, 0, MPI_COMM_WORLD);
}

void gather_spawned_organisms_other_ranks(environment* env, int rank, int size, MPI_Datatype organism_type) {
    // TODO: every process has its own list of spawned organisms, they need to be merged into a
    // single list in one process before environment_process_spawned_organisms() is called from
    // that process. Then, the list of spawned organisms should be sent back to the other
    // processes so organisms can be spawned in the correct locations. At the end, each process
    // should clear its own list of spawned organisms after the merge by doing
    // env.spawned[i].size = 0; This code should probably be placed in a function.

    // spawn new organisms
    // TODO: need to move all spawned organisms to the primary process before calling environment_process_spawned_organisms()
    environment_combine_spawned_organisms(env); // this sums over all threads in the current process, not across processes
    MPI_Gather(&env->spawned[0].size, 1, MPI_INT, NULL, 1, MPI_INT, 0, MPI_COMM_WORLD);

    spawned_org_array* spawned = &env->spawned[0];

    MPI_Gatherv(spawned->data, spawned->size, organism_type,
        spawned->data, NULL, NULL, organism_type,
        0, MPI_COMM_WORLD);
    
    // TODO: Scatter sizes and spawned organisms to all processes
    MPI_Scatter(NULL, 1, MPI_INT, &env->spawned[0].size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    spawned_org_array_ensure_capacity(spawned, env->spawned[0].size);

    MPI_Scatterv(env->spawned[0].data, NULL, NULL, organism_type, env->spawned[0].data, env->spawned[0].size, organism_type, 0, MPI_COMM_WORLD);
}


int main(int argc, char* argv[]) {
    // TODO: initialize MPI, get rank & size
    int provided;
    MPI_Init_thread(&argc, &argv,
        MPI_THREAD_FUNNELED, &provided);
    if (provided != MPI_THREAD_FUNNELED){
        printf("Error: MPI init failed\n");
        MPI_Finalize();
        return 0;
    }

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); // retrieve the current process' rank
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size); // retrieve the total number of processes

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
            fprintf(stderr, "usage: %s [-n num-iterations] [-r random-seed] [-m mutation-prob] [-d] [-t threads_per_proc] input output\n", argv[0]);
            return 1;
        }
    }
    if (optind + 2 != argc || iterations == 0 || num_threads == 0) {
        fprintf(stderr, "usage: %s [-n num-iterations] [-r random-seed] [-m mutation-prob] [-d] [-t threads_per_proc] input output\n", argv[0]);
        return 1;
    }
    const char* input_path = argv[optind];
    const char* output_path = argv[optind+1];

    // load the environment
    // NOTE: all processes load the same initial set of data
    environment env;
    if (!environment_load_from_path(&env, input_path, random_seed, num_threads)) { return 1; }
    env.mutation_prob = mutation_prob;
    size_t world_size = env.world_size;

    if (organism_debug) {
        organism_print(&env.grid[(world_size*world_size+world_size)/2]); // prints middle grid cell
    }

    // TODO: initialize the spawned organism MPI datatype
    MPI_Datatype organism_type;
    MPI_Type_create_struct(3, organism_field_counts, organism_field_offsets, organism_field_types, &organism_type);
    MPI_Type_commit(&organism_type);

    // TODO: have every process figure out which rows/cells it is responsible for
    int start_row = rank * world_size / size;
    int end_row = (rank + 1) * world_size / size;

    
    
    // start the timer
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // run the simulation
    for (size_t t = 0; t < iterations; t++) {
        // TODO: run the inner loop in parallel, but only over ones for this MPI process
        #pragma omp parallel for default(none) num_threads(num_threads) firstprivate(world_size, start_row, end_row) shared(env) schedule(dynamic, 16)
        for (size_t i = start_row * world_size; i < end_row * world_size; i++) {
            organism_update(&env.grid[i]);
        }
        
        // TODO: need to move all spawned organisms to the primary process before calling environment_process_spawned_organisms()
        if (rank == 0) {
            gather_spawned_organisms_rank_zero(&env, size, organism_type);
        } else {
            gather_spawned_organisms_other_ranks(&env, rank, size, organism_type);
        }

        // spawn new organisms
        for (size_t i = 0; i < env.spawned[0].size; i++) {
            // TODO: need to spawn the organisms in the correct processes
            organism_spawn(&env, &env.spawned[0].data[i]);
        }
        // TODO: clear the spawned organisms in all non-primary processes this needs to set the size to 0 for all spawned arrays
        env.spawned[0].size = 0;

        if (organism_debug && t % (1024*1024) == 0) { environment_print(&env); }
    }

    // TODO: sum number of tasks from all processes for distributed version
    environment_combine_tasks(&env); // this sums over all threads in the current process, not across processes

    if (rank != 0) {
        // TODO: move organisms to the primary process
    } else {
        // TODO: receive all of the organisms from other processes and place in our grid

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
    }

    // cleanup
    environment_free(&env);

    // TODO: cleanup/finalize anything added for MPI
    MPI_Type_free(&organism_type);
    MPI_Finalize();
    return 0;
}
