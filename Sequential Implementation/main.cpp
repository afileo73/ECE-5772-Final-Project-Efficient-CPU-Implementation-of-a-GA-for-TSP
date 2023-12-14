#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "consts.cpp"

#ifdef EMBEDDED
  #include <sys/time.h>
  #include <pthread.h>
#else
  #include <time.h>
#endif

#ifdef PARALLEL
  #include "GA_functions_parallel.cpp"
#else
  #include "GA_functions.cpp"
#endif

#ifdef PARALLEL
// Structure for thread arguments
typedef struct {
  int **pop;
  float *cost;
  int *parents;
  float **cost_table;
  int start;
  int end;
} TH_args;
#endif

int main(int argc, char **argv){
  // -------------Initialization-------------

  #ifdef TIMING
  #ifdef EMBEDDED
  // Timing Variables:
  struct timeval start, end, startgen, endgen;
  long t_us;
  gettimeofday(&start, NULL);
  #else
  float milliseconds;
  clock_t start, end, startgen, endgen;

  start = clock();
  #endif
  #endif

  #ifdef PARALLEL
    pthread_t *thread;
    TH_args *thread_args;
    thread_args = (TH_args *)calloc(NUM_THREADS, sizeof(TH_args));
    thread = (pthread_t *) malloc(NUM_THREADS*sizeof(pthread_t))
  #endif

  // Variable Initialization:
  int **pop; // The population
  // Allocate memory for each member of the populations chromosome
  pop = (int**) calloc(POPULATION_SIZE, sizeof(int*));
  int i, j;
  for(i = 0; i<POPULATION_SIZE; i++){
    // Allocate memory for each chromosome's genes
    pop[i] = (int*) calloc(NUM_CITIES, sizeof(int));
  }
  float *cost; // Each chromosomes cost
  cost = (float*)calloc(POPULATION_SIZE, sizeof(float));
  int *parents; // Selected parents to create next generation
  parents = (int*)calloc(POPULATION_SIZE*2, sizeof(int));
  float **cost_table;
  cost_table = (float**)calloc(NUM_CITIES,sizeof(float*));
  for(i = 0; i<NUM_CITIES; i++){
    cost_table[i] = (float*)calloc(NUM_CITIES,sizeof(float));
  }

  #ifdef PARALLEL
    // The range of the population a single thread should handle, rounded up
    int thread_range = (int)(((float)POPULATION_SIZE / (float)NUM_THREADS) + 0.5); 
    // Initializing thread arguments
    for(i=0; i < NUM_THREADS; i++){
      thread_args[i].pop = pop;
      thread_args[i].cost = cost;
      thread_args[i].parent = parents;
      thread_args[i].cost_table = cost_table;
      thread_args[i].start = (i * thread_range);
      thread_args[i].end = ((i+1)* thread_range) - 1;
      if (thread_args[i].end >= POPULATION_SIZE){
        thread_args[i].end = POPULATION_SIZE-1;
      }
    }
  #endif

  // Build Cost Table
  build_cost_table(cost_table);

  // Initialize Population
  initialize_population(pop);
  #ifdef DEBUG
    for(i=0; i<POPULATION_SIZE; i++){
      for(j=0; j<NUM_CITIES; j++){
        if(pop[i][j] < 0 || pop[i][j] >= NUM_CITIES){
          printf("Population member %d has invalid city %d at index %d\n", i, pop[i][j], j);
        }
      }
    }
  #endif

  // Cost Evaluation
  cost_update(pop, cost, cost_table);
  double min_cost = findleastcost(cost, cost_table);

  #ifdef TIMING
    #ifdef EMBEDDED
      gettimeofday(&end, NULL);
      t_us = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec;
      printf("Initialization took %ld us\n", tus);
    #else
      end = clock();
      milliseconds = ((double)(end - start)) / CLOCKS_PER_SEC;
      printf("Initialization took %f ms\n", milliseconds);
    #endif
  #endif

  #ifdef VERBOSE
    printf("Initial population least cost: %.0f\n", min_cost);
  #endif
  // -----------End Initialization-----------
  // -------------Begin GA Loop--------------
  bool stopping_criteria_met = false;
  int generation_count = 1;

  while(!stopping_criteria_met){
    #ifdef TIMING
      #ifdef EMBEDDED
      gettimeofday(&startgen, NULL);
        gettimeofday(&start, NULL);
      #else
        startgen = clock();
        start = clock();
      #endif
    #endif
    // Select Parents
    selection(cost, parents);

    #ifdef TIMING
      #ifdef EMBEDDED
        gettimeofday(&end, NULL);
        t_us = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec;
        printf("Gen %d: Selection took %ld us\n", generation_count,tus);
      #else
        end = clock();
        milliseconds = ((double)(end - start)) / CLOCKS_PER_SEC;
        printf("Gen %d: Selection took %f ms\n",generation_count, milliseconds);
      #endif
    #endif
    #ifdef DEBUG
      for(i=0; i<POPULATION_SIZE*2; i++){
        if(parents[i] < 0 || parents[i] >= POPULATION_SIZE){
          printf("Parent %d has invalid value %d\n", i, parents[i]);
        }
      }
    #endif
    #ifdef TIMING
      #ifdef EMBEDDED
        gettimeofday(&start, NULL);
      #else
        start = clock();
      #endif
    #endif
    // Crossover
    crossover(pop, parents, cost_table);

    #ifdef TIMING
      #ifdef EMBEDDED
        gettimeofday(&end, NULL);
        t_us = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec;
        printf("Gen %d: Crossover took %ld us\n",generation_count, tus);
      #else
        end = clock();
        milliseconds = ((double)(end - start)) / CLOCKS_PER_SEC;
        printf("Gen %d: Crossover took %f ms\n",generation_count, milliseconds);
      #endif
    #endif
    #ifdef DEBUG
      for(i=0; i<POPULATION_SIZE; i++){
        for(j=0; j<NUM_CITIES; j++){
          if(pop[i][j] < 0 || pop[i][j] >= NUM_CITIES){
            printf("Population member %d has invalid city %d at index %d\n", i, pop[i][j], j);
          }
        }
      }
    #endif

    #ifdef TIMING
      #ifdef EMBEDDED
        gettimeofday(&start, NULL);
      #else
        start = clock();
      #endif
    #endif
    // Mutation
    mutation(pop);
    #ifdef TIMING
      #ifdef EMBEDDED
        gettimeofday(&end, NULL);
        t_us = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec;
        printf("Gen %d: Mutation took %ld us\n",generation_count, tus);
      #else
        end = clock();
        milliseconds = ((double)(end - start)) / CLOCKS_PER_SEC;
        printf("Gen %d: Mutation took %f ms\n",generation_count, milliseconds);
      #endif
    #endif

    #ifdef DEBUG
      for(i=0; i<POPULATION_SIZE; i++){
        for(j=0; j<NUM_CITIES; j++){
          if(pop[i][j] < 0 || pop[i][j] >= NUM_CITIES){
            printf("Population member %d has invalid city %d at index %d\n", i, pop[i][j], j);
          }
        }
      }
    #endif

    // Population's Fitness
    #ifdef TIMING
      #ifdef EMBEDDED
        gettimeofday(&start, NULL);
      #else
        start = clock();
      #endif
    #endif
    cost_update(pop, cost, cost_table);
    min_cost = findleastcost(cost, cost_table);
    #ifdef TIMING
      #ifdef EMBEDDED
        gettimeofday(&end, NULL);
        t_us = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec;
        printf("Gen %d: Fitness took %ld us\n", generation_count,tus);
      #else
        end = clock();
        milliseconds = ((double)(end - start)) / CLOCKS_PER_SEC;
        printf("Gen %d: Fitness took %f ms\n",generation_count, milliseconds);
      #endif
    #endif
    printf("Generation %d's minimum cost: \t %.0f\n",generation_count,min_cost);
    // Stopping Conditions
    if(generation_count == NUM_GENERATIONS){
      stopping_criteria_met = true;
    }
    generation_count++;
    #ifdef TIMING
      #ifdef EMBEDDED
        gettimeofday(&endgen, NULL);
        t_us = (endgen.tv_sec - startgen.tv_sec)*1000000 + end.tv_usec - start.tv_usec;
        printf("Gen %d took %ld us\n",generation_count, tus);
      #else
        endgen = clock();
        milliseconds = ((double)(endgen - startgen)) / CLOCKS_PER_SEC;
        printf("Gen %d took %f ms\n",generation_count, milliseconds);
      #endif
    #endif
  }
  // --------------End GA Loop---------------
  // Free memory
  for(i=0; i<POPULATION_SIZE; i++){
    free(pop[i]);
  }
  free(pop);
  free(cost);
  free(parents);
  for(i = 0; i<NUM_CITIES; i++){
    free(cost_table[i]);
  }
  free(cost_table);
  #ifdef PARALLEL
    free(thread_args);
    free(thread);
  #endif

  return 0;
}