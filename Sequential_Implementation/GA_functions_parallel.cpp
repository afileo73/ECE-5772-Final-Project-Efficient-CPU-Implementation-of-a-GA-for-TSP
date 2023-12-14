#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <limits.h>
#include <math.h>
#include "consts.cpp"
#ifdef EMBEDDED
  #include <pthread.h>
#endif
#pragma once

// Structure for thread arguments
typedef struct {
  int **pop;
  float *cost;
  int *parents;
  float **cost_table;
  int start;
  int end;
  int seed;
} TH_args;

// Finds the linear distance between 2D coordinates
float L2distance(float x1, float y1, float x2, float y2) {
    float x_d = pow(x1 - x2, 2);
    float y_d = pow(y1 - y2, 2);
    return sqrt(x_d + y_d);
}

// Initialization cost table from the (x,y) locations of the cities
// Could be optimizied since distances are bidirectional, ie cost_table[i][j] = cost_table[j][i] 
void build_cost_table(float **cost_table){
  int k, j;
  for (k = 0; k < NUM_CITIES; k++) {
      for (j = 0; j < NUM_CITIES; j++) {
          if (k != j) {
              cost_table[k][j] = L2distance(city_x[k], city_y[k], city_x[j], city_y[j]);
          }
          else {
              cost_table[k][j] = 0.0;
          }
      }
  }
};

// Function for initializing
// each member of the population with a random permutation of the cities
void initialize_population(int **pop){
  int i, j;
  for(i = 0; i < POPULATION_SIZE; i++){
    // Give each gene a default value equal to its position in the array
    for(j = 0; j<NUM_CITIES; j++){
      pop[i][j] = j;
    }

    // Randomly shuffle each chromosomes gene positions
    // Skip index 0 since the first city is always the same
    int temp, pos;
    for(j = 1; j<NUM_CITIES; j++){
      pos = rand() % NUM_CITIES;
      temp = pop[i][j];
      pop[i][j] = pop[i][pos];
      pop[i][pos] = temp;
    }
  }
  return;
}

// Updates the cost of all chromosomes
void cost_update(int **pop, float *cost, float** cost_table){
  int i, j;

  // Evaluate every member of the population
  for(i = 0; i<POPULATION_SIZE; i++){
    cost[i] = 0.0; // Base cost
    // Loop through current chromosome and total cost
    for(j = 1; j<NUM_CITIES; j++){
      cost[i] += cost_table[pop[i][j-1]][pop[i][j]];
    }
  }
}

// Find the fittest member of the population
double findleastcost(float *cost, float** cost_table){
  int i;
  float minimum = cost[0];
  for(i = 1; i<POPULATION_SIZE; i++){
    if (cost[i] < minimum){
      minimum = cost[i];
    }
  }
  return(minimum);
};

// Perform a series of tournament selections to choose parents for the next
// generation of solutions
void* selection(void *slice){
  TH_args args = *( (TH_args *) slice); // 'slice' is a pointer to a structure

  float *cost = args.cost;
  int *parents = args.parents;
  int start = args.start;
  int end = args.end;
  start *= 2;
  end *= 2;
  //printf("Selection thread starts at %d ends at %d", start, end);
  int i, j;
  int * tournament;
  int temp_index, best_index;
  tournament = (int*)calloc(TOURNAMENT_SIZE, sizeof(int));
  // Select a two parents for every member of the next generation
  for(i = start; i != end; i++){
    // Select TOURNAMENT_SIZE number of the population to compete
    // in the tournament
    best_index = rand_r(&args.seed) % POPULATION_SIZE;
    for(j = 1; j < TOURNAMENT_SIZE; j++){
      temp_index = rand_r(&args.seed) % POPULATION_SIZE;
      if(cost[temp_index] < cost[best_index]){
        best_index = temp_index;
      }
    }
    parents[i] = best_index;
  }
}

// Subfunction of crossover for next valid index of parents
int getValidNextCity(int *parent, int *child, int current_index, bool* used_cities){
  int i, j;

  // Find the next valid city of the parent from the starting index
  for(i = current_index; i < NUM_CITIES; i++){
    // Check if city has been used in child yer
    if(!used_cities[parent[i]]){
      return(parent[i]);
    }
  }
  
  // If no valid city was found, find sequentially the next valid city in the child
  for(i = 0; i < NUM_CITIES; i++){
    if(!used_cities[i]){
      return(i);
    }
  }
  printf("You shouldn't be here! No valid city found in getValidNextCity()");
  return -1;
}

// Combine parents into children
void* crossover(void *slice){
  TH_args args = *((TH_args *) slice);
  int** pop = args.pop;
  int* parents = args.parents;
  float** cost_table = args.cost_table;
  int start = args.start;
  int end = args.end;

  int **new_pop; // The population
  // Allocate memory for each member of the populations chromosome
  new_pop = (int**) calloc(POPULATION_SIZE, sizeof(int*));
  int i, j;
  for(i = start; i!=end; i++){
    // Allocate memory for each chromosome's genes
    new_pop[i] = (int*) calloc(NUM_CITIES, sizeof(int));
  }
  // Produce a new child to replace every member of the populations
  for(i=start; i!=end; i++){
    new_pop[i][0] = 0; //First city is always zero
    bool used_cities[NUM_CITIES] = {};
    used_cities[0] = true;
    for(j = 1; j < NUM_CITIES; j++){
      int choice1 = getValidNextCity(pop[parents[i]],new_pop[i],j,used_cities);
      int choice2 = getValidNextCity(pop[parents[i+POPULATION_SIZE]],new_pop[i],j,used_cities);
      // Pick the better choice based on cost
      if(cost_table[new_pop[i][j-1]][choice1] < cost_table[new_pop[i][j-1]][choice2]){
        new_pop[i][j] = choice1;
        used_cities[choice1] = true;
      }else{
        new_pop[i][j] = choice2;
        used_cities[choice2] = true;
      }
    }
  }

  // Copy New Pop Data
  for(i=start; i!=end; i++){
    for(j=0; j<NUM_CITIES; j++){
      pop[i][j] = new_pop[i][j];
    }
  }

  // Free Memory
  for(i = start; i!=end; i++){
    free(new_pop[i]);
  }
  free(new_pop);
}

// Mutate random members of the population
void* mutation(void *slice){
  TH_args args = *( (TH_args *) slice); // 'slice' is a pointer to a structure
  int **pop = args.pop;
  int start = args.start;
  int end = args.end;

  int i;
  for(i = start; i != end; i++){
    // If a random percent chance occurs
    if((rand_r(&args.seed) % 100) <= MUTATION_CHANCE){
      // Select two random indexs and swap
      int index1 = rand() % NUM_CITIES;
      int index2 = rand() % NUM_CITIES;
      int temp = pop[i][index1];
      pop[i][index1] = pop[i][index2];
      pop[i][index2] = temp;
    }
  }
}