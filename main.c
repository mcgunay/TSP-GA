#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

const int BIGINT = 9999999;
enum CrossOverType {OX = 0, SCX = 1};
enum MutationType {ISM = 0, IVM = 1, SM = 2, Random = 3};

typedef struct node {
    int val;
    struct node * next;
} node_t;

struct offsprings{
    struct city* offspring1;
    struct city* offspring2;
};

struct city{
    int city_id;
    float x_coordinate;
    float y_coordinate;
};

struct individual{
    struct city* cities;
    int city_count;
};

struct population{
    struct individual* individuals;
    int population_count;
};

float calculate_distance(struct city* city1, struct city* city2){

    int x_distance = abs(city1->x_coordinate - city2->x_coordinate);
    int y_distance = abs(city1->y_coordinate - city2->y_coordinate);
    float distance = sqrt( (x_distance*x_distance) + (y_distance*y_distance) );

    return distance;
}

float** create_distance_matrix(struct city* cities, int dimension){
    float** distance_m = malloc(dimension * sizeof(float*));

    for (int i = 0; i < dimension; i++) {
        distance_m[i] = malloc(dimension * sizeof(float));

        for (int j = 0; j < dimension; j++) {

            if (i == j) {
                distance_m[i][j] = BIGINT;
                continue;
            }

            //float tmp = (float) (pow((double) cities[i][0] - c[j][0], 2.0) + pow((double) c[i][1] - c[j][1], 2.0));
            distance_m[i][j] = calculate_distance(&cities[i], &cities[j]);
        }
    }

    return distance_m;

}

void initialize_with_nearest_neighbor(struct population* pop, struct city* cities,float** distance_m, double percentage){

    int create_n_random = pop->population_count * percentage;
    srand(time(NULL));

    for(int i = 25; i< 25 + create_n_random; ++i){
        int random = rand() % 280;
        //Put first random city in the list
        pop->individuals[i].cities[0] = cities[random];

        node_t * head = NULL;
        head = malloc(sizeof(node_t));
        head->val = random;
        head->next = NULL;
        node_t* last = head;
        for(int j = 1; j< 280 ;j++){
            int nearest_neighbor_distance= BIGINT;
            int nearest_neighbor=0;
            for(int x = 0; x < 280; ++x){
                int picked = 0;
                if(distance_m[random][x] < nearest_neighbor_distance){
                    node_t * current = head;
                    while (current != NULL) {
                        if(x == current->val){
                           picked = 1;
                        }
                        current = current->next;
                    }
                    if(0 == picked){
                        nearest_neighbor_distance = distance_m[random][x];
                        nearest_neighbor = x;
                    }

                }
            }
            pop->individuals[i].cities[j] = cities[nearest_neighbor];
            random = nearest_neighbor;
            node_t* next = malloc(sizeof(node_t));
            next->val = nearest_neighbor;
            last->next = next;
            next->next = NULL;
            last = next;

        }
    }
}

void initialize_randomly(struct population* pop, struct city* cities, double percentage){

    int create_n_random = pop->population_count * percentage;
    srand(time(NULL));

    for(int i = 0; i< create_n_random; i++){
        //int random = rand() % 280;
        int* r = malloc(280 * sizeof(int));

        for(int i=0;i<280;++i){
            r[i]=i;
        }

        for (int i = 280-1; i >= 0; --i){
            //generate a random number [0, n-1]
            int j = rand() % (i + 1);

            //swap the last element with element at random index
            int temp = r[i];
            r[i] = r[j];
            r[j] = temp;
        }

        for(int j = 0; j< 280 ;j++){
            pop->individuals[i].cities[j] = cities[r[j]];

        }
        free(r);
    }
}

float calculate_fitness(struct city* cities, float** distance_m, int dimension){
    float fitness_value = 0;
    for(int i = 0;i<dimension - 1;i++){
        fitness_value = fitness_value + distance_m[cities[i].city_id][cities[i+1].city_id];
       // printf("fitness value: %f", fitness_value);
    }
    fitness_value = fitness_value + distance_m[cities[279].city_id][cities[0].city_id];

    return fitness_value;
}

struct offsprings* perform_order_crossover(struct city* parent1, struct city* parent2, int length, int dimension){

    struct city* child1 = malloc(dimension* sizeof(struct city));
    struct city* child2 = malloc(dimension* sizeof(struct city));

    srand(time(NULL));
//    int crossover_start = rand() % (dimension/2);
//    int crossover_end = crossover_start + dimension / 2;
    int crossover_start = 2;
    int crossover_end = 6;

    //First copy subset from parents to children
    for(int i = crossover_start;i<crossover_end;++i){
        child1[i] = parent2[i];
        child2[i] = parent1[i];
    }

    for(int i = crossover_end; i< dimension; ++i){
        int exists_in_subset = 0;
        for(int i1 = crossover_start; i1 < crossover_end; ++i1){
            if(parent2[i].city_id == child2[i1].city_id)
                exists_in_subset = 1;
        }
        if(0 == exists_in_subset)
            child2[i] = parent2[i];

        exists_in_subset = 0;
        for(int i1 = crossover_start; i1 < crossover_end; ++i1){
            if(parent1[i].city_id == child1[i1].city_id)
                exists_in_subset = 1;
        }
        if(0 == exists_in_subset)
            child1[i] = parent1[i];

    }

    int parent1_index = 0;
    int parent2_index = 0;

    for(int i = 0; i<dimension;++i){

        for(int j = parent2_index;j<dimension;++j){
            int exists_in_subset = 0;
            for(int i1 = 0; i1 < dimension; ++i1){
                if(parent2[j].city_id == child2[i1].city_id)
                    exists_in_subset = 1;
            }
            if(0 == exists_in_subset && child2[i].city_id == 0){
                child2[i] = parent2[j];
                parent2_index = j;
                break;
            }
        }
    }


    for(int i = 0; i<dimension;++i){

        for(int j = parent1_index;j<dimension;++j){
            int exists_in_subset = 0;
            for(int i1 = 0; i1 < dimension; ++i1){
                if(parent1[j].city_id == child1[i1].city_id)
                    exists_in_subset = 1;
            }
            if(0 == exists_in_subset && child1[i].city_id == 0){
                child1[i] = parent1[j];
                parent1_index = j;
                break;
            }
        }
    }

    struct offsprings* offsprings = malloc(sizeof(struct offsprings));
    offsprings->offspring1 = malloc(280* sizeof(struct city));
    offsprings->offspring2 = malloc(280* sizeof(struct city));

    offsprings->offspring1 = child1;
    offsprings->offspring2 = child2;

    return offsprings;
}

int main(int argc, char *argv[]){
    //Parameters are set via arg
    int iteration_count = 0;
    enum CrossOverType xover_type;
    enum MutationType mutation_type;

    //variables to read from file
    FILE *fp;
    char line[200];
    struct city* cities = (struct city*)malloc( 280 * sizeof(struct city));
    fp = fopen("../a280.tsp", "r"); // read mode

    if (fp == NULL)
    {
        perror("Error while opening the file.\n");
        exit(EXIT_FAILURE);
    }

    //Read arguments for algorithm
    if(argc == 4){
        if(argv[1] != NULL){
            sscanf(argv[1], "%d", &iteration_count);
            printf("iteration count = %d", iteration_count);
        }
        if(argv[2] != NULL){
            sscanf(argv[2], "%d", &xover_type);
        }
        if(argv[3] != NULL){
            sscanf(argv[3], "%d", &mutation_type);
        }
    }

    //Distance Matrix
    float** distance_m;

    //Pass cities to struct array
    int i = 0;
    int pass_first_lines = 0;
    int dimension = 0;
    int dimension_m = 0;
    const char* section = (char*)malloc(20* sizeof(char));
    char* s_dimension = "DIMENSION:";
    while ( ( fgets ( line, sizeof ( line), fp)) != NULL) {
        if(pass_first_lines++ < 6){
            if(sscanf ( line, "%s %d"
                    , section, &dimension) == 2) {
                if(strcmp(section,s_dimension) == 0){
                    dimension_m = dimension;
                }
            }
            continue;
        }

        struct city *n= (struct city*)malloc(sizeof(struct city));

        if ( ( sscanf ( line, "%d %f %f"
                , &n->city_id, &n->x_coordinate, &n->y_coordinate)) == 3){
            cities[i] = *n;
            i++;
        }
    }
    //close file
    fclose(fp);
    printf("dimension: %d", dimension_m);

    //Fill Distance Matrix
    distance_m = create_distance_matrix(cities, dimension);

    //printf("dimension between 1 and 2 is: %f", distance_m[0][1]);
    //Initialize a population;
    struct population* pop = malloc(sizeof(struct population));
    pop->population_count = 50;
    pop->individuals = malloc(50 * sizeof(struct individual));
    for(int i = 0; i < pop->population_count; i++){
        pop->individuals[i].cities = malloc(280 * sizeof(struct city));
    }

    //Initialize half of the population randomly
    initialize_randomly(pop, cities, 0.5);
    //struct city* tour = (struct city*)malloc( 280 * sizeof(struct city));
    initialize_with_nearest_neighbor(pop, cities,distance_m, 0.5);

//    float fit = calculate_fitness(pop->individuals[30].cities, distance_m, 280);
//    printf("fitness value = %f", fit);

    struct city* cities1 = malloc(9 * sizeof(struct city));
    cities1[0].city_id = 1;
    cities1[1].city_id = 2;
    cities1[2].city_id = 3;
    cities1[3].city_id = 4;
    cities1[4].city_id = 5;
    cities1[5].city_id = 6;
    cities1[6].city_id = 7;
    cities1[7].city_id = 8;
    cities1[8].city_id = 9;

    struct city* cities2 = malloc(9 * sizeof(struct city));
    cities2[0].city_id = 5;
    cities2[1].city_id = 7;
    cities2[2].city_id = 4;
    cities2[3].city_id = 9;
    cities2[4].city_id = 1;
    cities2[5].city_id = 3;
    cities2[6].city_id = 6;
    cities2[7].city_id = 2;
    cities2[8].city_id = 8;

    struct offsprings* off = perform_order_crossover(cities1, cities2, 0, 9);

    for (int j = 0; j < 9; ++j) {
        printf("%d ", off->offspring1[j].city_id);
    }
    printf("\n");
    for (int j = 0; j < 9; ++j) {
        printf("%d ", off->offspring2[j].city_id);
    }


    return 0;


    printf("Hello, World!\n");
    return 0;


}