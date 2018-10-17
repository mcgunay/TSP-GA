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

struct next_legit_genes{
    int next_of_first_parent;
    int next_of_second_parent;
};

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

struct offsprings* perform_OX(struct city* parent1, struct city* parent2, int length, int dimension){

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

int is_next_legit(struct city* parent, struct city* child, int current_index, int dimension){

    if(current_index  == dimension -1)
        return 0;

    int is_next_legit = 1;
    int next = parent[current_index + 1].city_id;
    for (int i = 0; i < dimension; ++i) {
        if(next == child[i].city_id){
            is_next_legit = 0;
            break;
        }
    }

    return is_next_legit;
}
struct next_legit_genes* find_next_legit_genes(struct city* parent1, struct city* parent2,int p1_index, int p2_index,int dimension,int child_size, struct city* offspring){

    int index_in_first_parent_next = 0;
    int index_in_second_parent_next = 0;

    int is_next_legit_p1 = is_next_legit(parent1, offspring, p1_index, dimension);
    int is_next_legit_p2 = is_next_legit(parent2, offspring, p2_index, dimension);

    //Next node is none so find me legit one!
    int found = 0;
    if(0 == is_next_legit_p1){
        for (int i = 1; i <=dimension; ++i) {
            for(int i1 = 0;i1 < child_size + 1;++i1){
                if(i == offspring[i1].city_id){
                    found = 1;
                    break;
                }
            }
            if(0 == found){
                for (int j = 0; j < dimension; ++j) {
                    if(i == parent1[j].city_id){
                        index_in_first_parent_next = j;
                        index_in_second_parent_next  = p2_index + 1;
                        break;
                    }

                }
                break;
            }
            found = 0;
        }
    }else if(0 == is_next_legit_p2){
        for (int i = 1; i <= dimension; ++i) {
            for(int i1 = 0;i1 < child_size + 1;++i1){
                if(i == offspring[i1].city_id){
                    found = 1;
                    break;
                }
            }
            if(0 == found){
                for (int j = 0; j < dimension; ++j) {
                    if(i == parent2[j].city_id){
                        index_in_second_parent_next = j;
                        index_in_first_parent_next  = p1_index+ 1;
                        break;
                    }

                }
                break;
            }
            found = 0;
        }
    }else{
        index_in_first_parent_next  = p1_index + 1;
        index_in_second_parent_next  = p2_index + 1;
    }

    struct next_legit_genes* next = malloc(sizeof(struct next_legit_genes));
    next->next_of_first_parent = index_in_first_parent_next;
    next->next_of_second_parent = index_in_second_parent_next;

    return next;
}
struct offsprings* perform_SCX(struct city* parent1, struct city* parent2, float** distance_m, int dimension){

    int size = 0;
    int index = 0;
    int cost_to_pick_from_p2 = 0;
    int cost_to_pick_from_p1 = 0;
    struct city* offspring = malloc(dimension * sizeof(struct city));
    offspring[0] = parent1[0];
    int current_city_id = offspring[0].city_id;
    int last_added_from_parent = 1;
    int last_added_city_index = 0;

    while(size < dimension - 1){

        int index_in_first_parent = 0;
        int index_in_second_parent = 0;
        int index_in_first_parent_next = 0;
        int index_in_second_parent_next = 0;

        if(1 == last_added_from_parent){
            for (int i = 0; i < dimension; ++i) {
                if(parent2[i].city_id == current_city_id){
                    index_in_second_parent = i;
                    index_in_first_parent = last_added_city_index;
                    break;
                }
            }
        }else{
            for (int i = 0; i < dimension; ++i) {
                if(parent1[i].city_id == current_city_id){
                    index_in_first_parent = i;
                    index_in_second_parent = last_added_city_index;
                    break;
                }
            }
        }

        struct next_legit_genes* next_indexes = find_next_legit_genes(parent1, parent2, index_in_first_parent, index_in_second_parent, dimension, index, offspring);

        index_in_first_parent_next = next_indexes->next_of_first_parent;
        index_in_second_parent_next = next_indexes->next_of_second_parent;

        free(next_indexes);

        if(parent1[index_in_first_parent_next].city_id == parent2[index_in_second_parent_next].city_id){
            offspring[index + 1] = parent1[index_in_first_parent_next];
            last_added_from_parent = 1;
            last_added_city_index = index_in_first_parent_next;
        }else{
            cost_to_pick_from_p1 = distance_m[parent1[index_in_first_parent].city_id - 1][parent1[index_in_first_parent_next].city_id - 1];

            cost_to_pick_from_p2 = distance_m[parent2[index_in_second_parent].city_id - 1][parent2[index_in_second_parent_next].city_id - 1];

            if(cost_to_pick_from_p2 < cost_to_pick_from_p1){
                offspring[index + 1] = parent2[index_in_second_parent_next];
                last_added_from_parent = 2;
                last_added_city_index = index_in_second_parent_next;
            }else{
                offspring[index + 1] = parent1[index_in_first_parent_next];
                last_added_from_parent = 1;
                last_added_city_index = index_in_first_parent_next;
            }
        }

        current_city_id = offspring[index + 1].city_id;
        index++;
        size++;
    }

    return offspring;

}

void get_two_different_random_number(int* number1, int* number2, int dimension){

     srand(time(NULL));

    do {
        *number1 = rand() % (dimension);
        *number2 = rand() % (dimension);
    }while(*number1 == *number2);

}

void perform_inversion_mutation(struct city* child,int dimension){

    int* random1 = malloc(sizeof(int));
    int* random2 = malloc(sizeof(int));
    int first_random_gene = 0, second_random_gene = 0;

    get_two_different_random_number(random1, random2, dimension);

    if(random1 < random2){
        first_random_gene = *random1;
        second_random_gene = *random2;
    }else{
        first_random_gene = *random2;
        second_random_gene = *random1;
    }
    free(random1);
    free(random2);

    first_random_gene = 1;
    second_random_gene = 4;
    struct city* temp = malloc(sizeof(struct city));

    int head = first_random_gene, tail = second_random_gene;
    while(head < tail){

        *temp = child[head];
        child[head] = child[tail];
        child[tail] = *temp;

        head++;
        tail--;
    }

    free(temp);
}

void perform_swap_mutation(struct city* child,int dimension){

    int* random1 = malloc(sizeof(int));
    int* random2 = malloc(sizeof(int));
    int first_random_gene = 0, second_random_gene = 0;

    get_two_different_random_number(random1, random2, dimension);

    if(random1 < random2){
        first_random_gene = *random1;
        second_random_gene = *random2;
    }else{
        first_random_gene = *random2;
        second_random_gene = *random1;
    }
    free(random1);
    free(random2);

    struct city* temp = malloc(sizeof(struct city));
     *temp = child[first_random_gene];
    child[first_random_gene] = child[second_random_gene];
    child[second_random_gene] = *temp;

    free(temp);
}

void perform_Insert_mutation(struct city* child,int dimension){

//    srand(time(NULL));
//    int random1 = 0;
//    int random2 = 0;
//
//    do {
//         random1 = rand() % (dimension);
//         random2 = rand() % (dimension);
//    }while(random1 == random2);
    int* random1 = malloc(sizeof(int));
    int* random2 = malloc(sizeof(int));

    get_two_different_random_number(random1, random2, dimension);
    int first_random_gene = 0, second_random_gene = 0;

    if(random1 < random2){
        first_random_gene = *random1;
        second_random_gene = *random2;
    }else{
        first_random_gene = *random2;
        second_random_gene = *random1;
    }
    free(random1);
    free(random2);

    int intermediate_array_size = second_random_gene - first_random_gene - 1;

    struct city* intermediate = malloc(intermediate_array_size* sizeof(struct city));

    memcpy(intermediate, child + first_random_gene + 1, (intermediate_array_size * sizeof(struct city)));

    memmove(child + first_random_gene + 1, child + second_random_gene, sizeof(struct city));

    memcpy(child + first_random_gene + 2, intermediate, (intermediate_array_size * sizeof(struct city)));

    free(intermediate);

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

    struct city* cities1 = malloc(7 * sizeof(struct city));
    cities1[0].city_id = 1;
    cities1[1].city_id = 5;
    cities1[2].city_id = 7;
    cities1[3].city_id = 3;
    cities1[4].city_id = 6;
    cities1[5].city_id = 4;
    cities1[6].city_id = 2;


    struct city* cities2 = malloc(7 * sizeof(struct city));
    cities2[0].city_id = 1;
    cities2[1].city_id = 6;
    cities2[2].city_id = 2;
    cities2[3].city_id = 4;
    cities2[4].city_id = 3;
    cities2[5].city_id = 5;
    cities2[6].city_id = 7;

//    float** distance_mm = malloc(7 * sizeof(float*));
//
//    for (int i = 0; i < dimension; i++) {
//        distance_mm[i] = malloc(7 * sizeof(float));
//    }
//
//
//    //struct offsprings* off = perform_OX(cities1, cities2, 0, 9);
//    struct city* ofs = perform_SCX(cities1, cities2, distance_mm, 7);
//
//    for (int j = 0; j < 7; ++j) {
//        printf("%d ", ofs[j].city_id);
//    }

    struct city* cities3 = malloc(9 * sizeof(struct city));
    cities3[0].city_id = 1;
    cities3[1].city_id = 2;
    cities3[2].city_id = 3;
    cities3[3].city_id = 4;
    cities3[4].city_id = 5;
    cities3[5].city_id = 6;
    cities3[6].city_id = 7;
    cities3[7].city_id = 8;
    cities3[8].city_id = 9;

    perform_inversion_mutation(cities3, 9);

        for (int j = 0; j < 9; ++j) {
        printf(" insert mutation %d \n", cities3[j].city_id);
    }

    return 0;


    printf("Hello, World!\n");
    return 0;


}


//distance_mm[0][0] = 9999;
//distance_mm[0][1] = 75;
//distance_mm[0][2] = 99;
//distance_mm[0][3] = 9;
//distance_mm[0][4] = 35;
//distance_mm[0][5] = 63;
//distance_mm[0][6] = 8;
//
//distance_mm[1][0] = 51;
//distance_mm[1][1] = 9999;
//distance_mm[1][2] = 86;
//distance_mm[1][3] = 46;
//distance_mm[1][4] = 88;
//distance_mm[1][5] = 29;
//distance_mm[1][6] = 20;
//
//distance_mm[2][0] = 100;
//distance_mm[2][1] = 5;
//distance_mm[2][2] = 9999;
//distance_mm[2][3] = 16;
//distance_mm[2][4] = 28;
//distance_mm[2][5] = 35;
//distance_mm[2][6] = 28;
//
//distance_mm[3][0] = 20;
//distance_mm[3][1] = 45;
//distance_mm[3][2] = 11;
//distance_mm[3][3] = 9999;
//distance_mm[3][4] = 59;
//distance_mm[3][5] = 53;
//distance_mm[3][6] = 49;
//
//distance_mm[4][0] = 86;
//distance_mm[4][1] = 63;
//distance_mm[4][2] = 33;
//distance_mm[4][3] = 65;
//distance_mm[4][4] = 9999;
//distance_mm[4][5] = 76;
//distance_mm[4][6] = 72;
//
//distance_mm[5][0] = 36;
//distance_mm[5][1] = 53;
//distance_mm[5][2] = 89;
//distance_mm[5][3] = 31;
//distance_mm[5][4] = 21;
//distance_mm[5][5] = 9999;
//distance_mm[5][6] = 52;
//
//distance_mm[6][0] = 58;
//distance_mm[6][1] = 31;
//distance_mm[6][2] = 43;
//distance_mm[6][3] = 67;
//distance_mm[6][4] = 52;
//distance_mm[6][5] = 60;
//distance_mm[6][6] = 9999;