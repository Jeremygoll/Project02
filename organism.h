/**
 * Defines the interface for an organism and the environment. Two of the
 * functions are not defined in organism.c and must be implemented elsewhere.
 * This file should mostly not need to be modifed with the exception being the
 * environment struct. 
 */
#pragma once

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>


////////////////////////////
///// Global Constants /////
////////////////////////////

// Organism Properties and Constants
#define MAX_GENOME_LENGTH 256 // no organism will be created with a genome larger than this (official implementation has 2048, but our program is more restricted)
#define MAX_STACK_SIZE 8   // maximum size of the stacks (official implementation has 10)
#define NUM_STACKS 2       // number of stacks
#define MAX_TEMPLATE_LEN 8 // maximum length of templates (official implementation has 10), also applies to the number of recent copied instructions that should be remembered)

// Definitions of the registers and heads in the organism
#define AX 0
#define BX 1
#define CX 2
#define NUM_REGISTERS 3
#define IP 0
#define READ_HEAD 1
#define WRITE_HEAD 2
#define FLOW_HEAD 3
#define NUM_HEADS 4

// Tasks to be completed
#define TASK_NOT     0
#define TASK_NAND    1
#define TASK_AND     2
#define TASK_OR_NOT  3
#define TASK_OR      4
#define TASK_AND_NOT 5
#define TASK_NOR     6
#define TASK_XOR     7
#define TASK_EQUAL   8
#define NUM_TASKS    9

// Energy they granted by tasks
static const uint32_t ENERGY_PER_TASK[] = {1024, 1024, 2048, 2048, 4096, 4096, 8192, 8192, 16384}; // last one is enough to fully replicate plus some


//////////////////////
///// Structures /////
//////////////////////

typedef struct _environment environment;


/**
 * An individual organism within the environment.
 * Each organism requires ~188 bytes more than just the double-length genome (total 512 with the MAX_GENOME_LENGTH of 256)
 */
typedef struct _organism {
    // NOTE: do NOT change this structure

    // information about the organism directly related to the CPU processing
    uint8_t memory[2*MAX_GENOME_LENGTH]; // cyclic memory buffer (CPU instructions and flags)
    uint32_t size; // used size of the memory buffer
    union {
        uint32_t heads[NUM_HEADS]; // "heads": instruction pointer, read, write, and flow heads
        struct { uint32_t ip, read_head, write_head, flow_head; }; // give names to the values in the array (sometimes array is easier to use)
    };
    union {
        int32_t registers[NUM_REGISTERS]; // registers AX, BX, and CX
        struct { int32_t ax, bx, cx; }; // give names to the values in the array (sometimes array is easier to use)
    };
    int32_t stacks[NUM_STACKS][MAX_STACK_SIZE]; // 2 stacks
    uint32_t stack_sizes[NUM_STACKS]; // the number of values currently in the stacks
    uint8_t stack_active; // which stack is active
    //void *in_buffer, *out_buffer; // these are not explicit within this implementation

    // additional information required for some instructions
    unsigned int random_seed; // the organism's personal random seed
    bool has_allocated; // if h-alloc has been called and h-divide hasn't been
    uint32_t energy; // the amount of energy the organism has, it takes 64 units to execute one instruction, gain 1 per update
    uint32_t genome_size; // original size of the genome
    uint32_t instructions_executed; // number of instructions executed (it's 'age')
    uint8_t recent_copies[MAX_TEMPLATE_LEN], recent_copies_offset, recent_copies_size; // recently copied values (cyclic buffer)
    uint32_t num_inputs; // number of inputs this organism has obtained
    uint32_t tasks_completed[NUM_TASKS]; // number of times this organism has completed each task

    // link to the environment
    uint32_t location_i, location_j; // the location within the environment
    environment* env; // the environment this organism is in
} organism;

/**
 * The information required to spawn an organism so it can be spawned at a later time.
 */
typedef struct _spawned_organism {
    // NOTE: do NOT change this structure
    uint8_t genome[MAX_GENOME_LENGTH];
    uint32_t size, loc;
} spawned_organism;

/**
 * A dynamic array of spawned organisms.
 */
typedef struct _spawned_org_array {
    spawned_organism* data;
    size_t size;
    size_t capacity;
} spawned_org_array;

/**
 * An array with a count for each task.
 */
typedef struct { size_t counts[NUM_TASKS]; } task_counts;

/**
 * The environment containing all of the organisms.
 * Any organism with a size of 0 is non-existent.
 */
struct _environment {
    // NOTE: you may change this structure to save additional information as needed
    organism* grid; // 2D array of organisms; world_size-by-world_size
    size_t world_size; // the size of the grid
    unsigned int num_threads; // number of threads being used
    double mutation_prob; // probability of mutation
    task_counts* tasks_completed; // number of times any organism has completed each task, one list for each thread
    spawned_org_array* spawned; // list of organisms to be spawned at the end of the iteration, one list for each thread
};


///////////////////////////////////////////////////////
///// Organism Info for Creating Custom MPI Types /////
///////////////////////////////////////////////////////
// To enable this section, you have to include mpi.h before you include this file
#ifdef MPI_VERSION
// These do not include the final field (the pointer to the environment)
static const int organism_field_counts[] = {
    MAX_GENOME_LENGTH*2, 1, NUM_HEADS, NUM_REGISTERS, NUM_STACKS*MAX_STACK_SIZE, NUM_STACKS, 1,
    1, 1, 1, 1, 1, MAX_TEMPLATE_LEN, 1, 1, 1, NUM_TASKS, 1, 1
};
static const MPI_Aint organism_field_offsets[] = {
    offsetof(organism, memory), offsetof(organism, size), offsetof(organism, heads), offsetof(organism, registers),
    offsetof(organism, stacks), offsetof(organism, stack_sizes), offsetof(organism, stack_active),
    offsetof(organism, random_seed), offsetof(organism, has_allocated), offsetof(organism, energy), offsetof(organism, genome_size),
    offsetof(organism, instructions_executed),
    offsetof(organism, recent_copies), offsetof(organism, recent_copies_offset), offsetof(organism, recent_copies_size),
    offsetof(organism, num_inputs), offsetof(organism, tasks_completed),
    offsetof(organism, location_i), offsetof(organism, location_j)
};
static const MPI_Datatype organism_field_types[] = {
    MPI_BYTE, MPI_UNSIGNED, MPI_UNSIGNED, MPI_INT, MPI_INT, MPI_UNSIGNED, MPI_BYTE,
    MPI_UNSIGNED, MPI_BYTE, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_BYTE, MPI_BYTE, MPI_BYTE, MPI_UNSIGNED, MPI_UNSIGNED,
    MPI_UNSIGNED, MPI_UNSIGNED
};

// Simplier spawned_organism - much easier to transfer over MPI
static const int spawned_organism_field_counts[] = { MAX_GENOME_LENGTH, 1, 1 };
static const MPI_Aint spawned_organism_field_offsets[] = { offsetof(spawned_organism, genome), offsetof(spawned_organism, size), offsetof(spawned_organism, loc) };
static const MPI_Datatype spawned_organism_field_types[] = { MPI_BYTE, MPI_UNSIGNED, MPI_UNSIGNED };
#endif


///////////////////////////////////////////////
///// Core Environment/Organism Functions /////
///////////////////////////////////////////////

/**
 * Load the environment from the path to a given text file. See
 * environment_load_from_file() for details.
 */
bool environment_load_from_path(environment* env, const char* path, unsigned int random_seed, int num_threads);

/**
 * Load the environment from the given text file. The file must have a single
 * number on the first line for the world size and every line after that has
 * 2 numbers separated by whitespace for the location of an organism followed
 * by the genome for the organism (which is a sequence of lowercase letters).
 * This returns false if there are any issues reading the file or the data is
 * not a proper environment (it also prints out an error message).
 * 
 * The random seed is used to set the random seed for each organism and the
 * number of threads is the number of separate lists of spawned organisms and
 * tasks to be allocated.
 */
bool environment_load_from_file(environment* env, FILE* file, unsigned int random_seed, int num_threads);

/**
 * Save the environment to the path of the given text file. See
 * environment_save_to_file() for details.
 */
bool environment_save_to_path(environment* env, const char* path);

/**
 * Save the environment to the given text file. It will have the same format as
 * used by environment_load_from_file().
 */
bool environment_save_to_file(environment* env, FILE* file);

/**
 * Setup the grid for an environment. This allocates the memory, clears the
 * grid making all cells unoccupied / all organisms non-existent, and sets up a
 * few values so the organisms link back to the environment. Additionally, this
 * sets every organism's random seed based on the given random seed. The
 * env->grid field must be freed after this is called. Before calling this
 * function, make sure world_size is set properly in the environment.
 */
bool environment_alloc_grid(environment* env, unsigned int random_seed);

/**
 * Free the memory used by the environment. This includes the grid of organisms,
 * the list of organisms to be spawned, and the list of tasks. This does not
 * free the environment itself.
 */
void environment_free(environment* env);

/**
 * Print out the environment. This displays it as a grid with spaces for empty
 * cells, single dots for occupied cells with weak (or the only) organisms, and
 * then blocks of various brightness for cells with stronger organisms.
 */
void environment_print(const environment* env);

/**
 * Set task counts for the environemnt to 0. Does not effect organisms.
 */
void environment_clear_tasks(environment* env);

/**
 * Combines all task counts for all threads into the first thread and clears
 * the counts for all other threads. This must be called from only one thread.
 */
void environment_combine_tasks(environment* env);

/**
 * Adds an organism to the list of organisms to be spawned in the environment
 * at the end of the current iteration. The genome is copied so after this is
 * called it can be freed without affecting the organism. The location is the
 * combined i*world_size+j value.
 * 
 * The environment has a separate list for each thread so that this function
 * does not need to be thread-safe.
 */
void environment_add_to_be_spawned(environment* env, uint8_t* genome, uint32_t size, uint32_t loc);

/**
 * Initializes a spawned organism array. This must be called before using the
 * array.
 */
spawned_org_array* spawned_org_array_init(spawned_org_array* arr);

/**
 * Resizes the list of spawned organisms to have at least the given capacity.
 * This function should only be called from a single thread at a time. It is
 * not thread-safe with any other function that modifies this spawned list.
 */
bool spawned_org_array_ensure_capacity(spawned_org_array* arr, size_t capacity);

/**
 * Combine the list of spawned organisms. This merges the list of spawned
 * organisms from all threads into a single list. They will end up in the thread 
 * 0's spawned list but no further processing.
 * 
 * This must be called from thread 0.
 */
void environment_combine_spawned_organisms(environment* env);

/**
 * Process the list of spawned organisms. This merges the list of spawned
 * organisms from all threads into a single list and then has them determine
 * the new locations for the organisms. They will end up in the thread 0's
 * spawned list and be sorted by their new location.
 * 
 * This must be called from thread 0.
 */
void environment_process_spawned_organisms(environment* env);

/**
 * Check if a cell in the environment is occupied. Used like:
 *      cell_occupied(&env->grid[i*ws + j]);
 */
bool cell_occupied(organism* org);

/**
 * Converts a genome from an alphabetical sequence (e.g. "wzcagcccccccccc...")
 * to one usable by an organism.
 */
void convert_genome(const char* genome_alpha, uint8_t* genome, uint32_t size);

/**
 * Converts a genome from an alphabetical sequence (e.g. "wzcagcccccccccc...")
 * to one usable by an organism. The memory is copied and the copy is returned.
 * The returned memory must be freed.
 */
uint8_t* convert_genome_copy(const char* genome_alpha, uint32_t size);

/**
 * Converts a genome to an alphabetical sequence (e.g. "wzcagcccccccccc...").
 * The given buffer for genome_alpha must have room for the null terminator.
 */
void convert_genome_to_alpha(const uint8_t* genome, char* genome_alpha, uint32_t size);

/**
 * Converts a genome to an alphabetical sequence (e.g. "wzcagcccccccccc...").
 * to one usable by an organism. The memory is copied and the copy is returned.
 * The returned memory must be freed.
 */
char* convert_genome_to_alpha_copy(const uint8_t* genome, uint32_t size);

/**
 * Set an organism to have the given genome. The genome must be converted.
 * The genome is copied so after this is called it can be freed without
 * affecting the organism. If there is already an organism it is erased/cleared.
 */
void organism_set(organism* org, uint8_t* genome, uint32_t size);

/**
 * Spawn an organism in the environment. The location is guaranteed to be
 * within the environment (from 0,0 to WORLD_SIZE,WORLD_SIZE). The genome must
 * be immediately used or copied as after calling this function that memory may
 * be overwritten.
 */
void organism_spawn(environment* env, spawned_organism* org);

/**
 * Update an organism by having it execute an instruction (if it has enough
 * energy) and then giving it a small amount of energy. The organism may die
 * (if it has grown too old) and it may divide. If the organism is non-existent
 * (it is cleared/the cell is not occupied) this function does nothing but is
 * still safe to call.
 */
void organism_update(organism* o);

/**
 * Print out all of the details about an organism. This includes its genome (as
 * an alphebetical string), the locations of the instruction point and heads,
 * the values of the registers, the stacks, recent inputs, number of
 * instructions executed, and number of tasks completed. If the organism is
 * non-existent (it is cleared/the cell is not occupied) "Non-existent" is
 * printed.
 */
void organism_print(const organism* org);

/**
 * Reset an organism. It only keeps:
 *   - genome (however any additional flags are removed)
 *   - energy
 *   - number of tasks completed
 * Everything else is as if it had be cleared and set.
 */
void organism_reset(organism* org);

/**
 * Clear an organism causing it to become non-existent. If this is a cell
 * within an environment, that cell is now un-occupied.
 */
void organism_clear(organism* org);
