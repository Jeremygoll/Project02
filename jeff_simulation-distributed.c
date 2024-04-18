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
#include <string.h>
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


#define MAX_MPI_SIZE 16  // max number of processes for MPI

/**
 * Process all spawned organisms for an MPI system.
 * This function is to be called by rank 0.
 */
void mpi_process_spawned_organisms_rank0(environment* env, MPI_Datatype MPI_Spawned_Org_Type) {
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
    int counts[size]; // counts of spawned organisms from/to each process
    int displacements[size]; // displacements for Gatherv/Scatterv (cumsum of counts)
    static size_t world_size = 0;
    static size_t end_cells[MAX_MPI_SIZE]; // number of rows for each process
    if (env->world_size != world_size) {
        world_size = env->world_size;
        for (int r = 0; r < size; r++) { end_cells[r] = (r + 1) * world_size / size * world_size; }
    }

    // gather all spawned organisms into the primary process
    spawned_org_array* array = &env->spawned[0];
    MPI_Gather(
        &array->size, 1, MPI_INT,
        counts, 1, MPI_INT, 0, MPI_COMM_WORLD
    );
    size_t sum = 0;
    for (int r = 0; r < size; r++) { displacements[r] = sum; sum += counts[r]; }
    spawned_org_array_ensure_capacity(array, sum);
    array->size = sum;
    MPI_Gatherv(
        MPI_IN_PLACE, array->size, MPI_Spawned_Org_Type,
        array->data, counts, displacements, MPI_Spawned_Org_Type, 0, MPI_COMM_WORLD
    );

    // process all spawned organisms
    environment_process_spawned_organisms(env);

    // figure out which spawned organisms need to be sent to each process 
    displacements[0] = 0;
    int cur_rank = 0;
    size_t n = array->size;
    spawned_organism* data = array->data;
    for (int i = 0; i < n; i++) {
        while (data[i].loc >= end_cells[cur_rank]) {
            cur_rank++;
            displacements[cur_rank] = i;
        }
    }
    for (int r = 0; r < size-1; r++) { counts[r] = displacements[r+1] - displacements[r]; }
    counts[size-1] = n - displacements[size-1];

    // distribute organisms from the primary process that need to be spawned by other processes
    MPI_Scatter(
        counts, 1, MPI_INT,
        &array->size, 1, MPI_INT, 0, MPI_COMM_WORLD
    );
    MPI_Scatterv(
        array->data, counts, displacements, MPI_Spawned_Org_Type,
        MPI_IN_PLACE, counts[0], MPI_Spawned_Org_Type, 0, MPI_COMM_WORLD
    );
}

/** Process all spawned organisms for an MPI system. This is run be all non-rank-0 processes. */
void mpi_process_spawned_organisms_worker(environment* env, MPI_Datatype MPI_Spawned_Org_Type) {
    // gather all spawned organisms into the primary process
    spawned_org_array* array = &env->spawned[0];
    MPI_Gather(
        &array->size, 1, MPI_INT,
        NULL, 0, MPI_INT, 0, MPI_COMM_WORLD
    );
    MPI_Gatherv(
        array->data, array->size, MPI_Spawned_Org_Type,
        NULL, NULL, NULL, MPI_Spawned_Org_Type, 0, MPI_COMM_WORLD
    );

    // distribute organisms from the primary process that need to be spawned by this process
    MPI_Scatter(
        NULL, 0, MPI_INT,
        &array->size, 1, MPI_INT, 0, MPI_COMM_WORLD
    );
    spawned_org_array_ensure_capacity(array, array->size);
    MPI_Scatterv(
        NULL, NULL, NULL, MPI_Spawned_Org_Type,
        array->data, array->size, MPI_Spawned_Org_Type, 0, MPI_COMM_WORLD
    );
}

int main(int argc, char* argv[]) {
    // TODO: initialize MPI, get rank & size
    int provided, rank, size;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED) {
        fprintf(stderr, "MPI_Init_thread did not provide the requested level of thread support\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    if (optind + 2 != argc || iterations == 0 || num_threads == 0) {
        fprintf(stderr, "usage: %s [-n num-iterations] [-r random-seed] [-m mutation-prob] [-d] [-t threads_per_proc] input output\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
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
    MPI_Datatype MPI_Spawned_Org_Type;
    MPI_Type_create_struct(sizeof(spawned_organism_field_counts)/sizeof(spawned_organism_field_counts[0]),
        spawned_organism_field_counts, spawned_organism_field_offsets, spawned_organism_field_types, &MPI_Spawned_Org_Type);
    MPI_Type_commit(&MPI_Spawned_Org_Type);

    // initialize the organism MPI datatype
    // MPI_Datatype MPI_Organism_Type;
    // MPI_Type_create_struct(sizeof(organism_field_counts)/sizeof(organism_field_counts[0]),
    //    organism_field_counts, organism_field_offsets, organism_field_types, &MPI_Organism_Type);
    // MPI_Type_commit(&MPI_Organism_Type);

    // TODO: have every process figure out which rows/cells it is responsible for
    if (world_size < size) { if (rank == 0) { fprintf(stderr, "cannot run more processes than rows in the grid\n"); } MPI_Abort(MPI_COMM_WORLD, 1); }
    size_t start_row = rank * world_size / size;
    size_t end_row = (rank + 1) * world_size / size;
    //size_t num_rows = end_row - start_row;
    size_t start_cell = start_row * world_size;
    size_t end_cell = end_row * world_size;
    const size_t max_rows_per_process = (world_size + size - 1) / size;
    //int rank_above = (rank+1)%size, rank_below = (rank==0?size:rank)-1;

    // start the timer
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // run the simulation
    for (size_t t = 0; t < iterations; t++) {
        // TODO: run the inner loop in parallel, but only over ones for this MPI process
        #pragma omp parallel for default(none) num_threads(num_threads) \
            firstprivate(start_cell, end_cell) shared(env) schedule(dynamic, 16)
        for (size_t i = start_cell; i < end_cell; i++) {
            organism_update(&env.grid[i]);
        }

        // TODO: every process has its own list of spawned organisms, they need to be merged into a
        // single list in one process before environment_process_spawned_organisms() is called from
        // that process. Then, the list of spawned organisms should be sent back to the other
        // processes so organisms can be spawned in the correct locations. At the end, each process
        // should clear its own list of spawned organisms after the merge by doing
        // env.spawned[i].size = 0; This code should probably be placed in a function.

        // combine all spawned organisms from all threads in this process to a single array 
        environment_combine_spawned_organisms(&env);

        // process all spawned organisms
        if (rank == 0) { mpi_process_spawned_organisms_rank0(&env, MPI_Spawned_Org_Type); }
        else { mpi_process_spawned_organisms_worker(&env, MPI_Spawned_Org_Type); }

        // spawn new organisms
        for (size_t i = 0; i < env.spawned[0].size; i++) {
            organism_spawn(&env, &env.spawned[0].data[i]);
        }
        // clear the spawned organisms in all non-primary processes this needs to set the size to 0 for all spawned arrays
        env.spawned[0].size = 0;

        if (organism_debug && t % (1024*1024) == 0) { environment_print(&env); }
    }

    // TODO: sum number of tasks from all processes for distributed version
    environment_combine_tasks(&env); // this sums over all threads in the current process, not across processes
    MPI_Reduce(rank == 0 ? MPI_IN_PLACE : env.tasks_completed, env.tasks_completed, NUM_TASKS, MPI_SIZE_T, MPI_SUM, 0, MPI_COMM_WORLD);

    spawned_organism* orgs = (spawned_organism*)malloc(max_rows_per_process*world_size*sizeof(spawned_organism));
    if (rank != 0) {
        // TODO: move organisms to the primary process
        for (int i = start_cell; i < end_cell; i++) {
            memcpy(orgs[i-start_cell].genome, env.grid[i].memory, env.grid[i].genome_size);
            orgs[i-start_cell].size = env.grid[i].genome_size;
            orgs[i-start_cell].loc = i;
        }
        MPI_Send(orgs, end_cell-start_cell, MPI_Spawned_Org_Type, 0, 0, MPI_COMM_WORLD);

    } else {
        // TODO: receive all of the organisms from other processes and place in our grid
        for (int r = 1; r < size; r++) {
            MPI_Status status;
            int count;
            MPI_Recv(orgs, max_rows_per_process*world_size, MPI_Spawned_Org_Type, r, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_Spawned_Org_Type, &count);
            for (size_t i = 0; i < count; i++) {
                organism_set(&env.grid[orgs[i].loc], orgs[i].genome, orgs[i].size);
            }
        }

        // final debug print
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
    free(orgs);
    MPI_Type_free(&MPI_Spawned_Org_Type);
    //MPI_Type_free(&MPI_Organism_Type);
    MPI_Finalize();
    return 0;
}
