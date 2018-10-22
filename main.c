#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sodium.h>

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

uint32_t get_one_random_number(int dimension){

    char myString[32];
    uint32_t myInt;
    if (sodium_init() < 0) {
        /* panic! the library couldn't be initialized, it is not safe to use */
        return 1;
    }
    /* myString will be an array of 32 random bytes, not null-terminated */
    randombytes_buf(myString, 32);
    //const unsigned char seed[randombytes_SEEDBYTES] = {"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"};
    //randombytes_buf_deterministic(myString, randombytes_SEEDBYTES, seed);

    /* myInt will be a random number between 0 and 9 */
    return randombytes_uniform(dimension);
}

void get_two_different_random_number(int* number1, int* number2, int dimension){

    //srand(time(NULL));

    do {
        *number1 = get_one_random_number(dimension);
        *number2 = get_one_random_number(dimension);
    }while(*number1 == *number2);

}

void get_n_different_random_number(int* numbers, size_t size, int dimension){

    //srand(time(NULL));
    int i = 0;
    int number = 0;

    while(i < size){
        int exist = 0;
        do{
            exist = 0;
            number = get_one_random_number(dimension);
            for (int j = 0; j < i; ++j) {
                if(number == numbers[j])
                    exist = 1;
            }
        }while(exist == 1);
        exist = 0;
        numbers[i] = number;
        i++;
    }

}

float calculate_distance(struct city* city1, struct city* city2){

    int x_distance = abs(city1->x_coordinate - city2->x_coordinate);
    int y_distance = abs(city1->y_coordinate - city2->y_coordinate);
    float distance = sqrt( (x_distance*x_distance) + (y_distance*y_distance) );

    return distance;
}

double calculate_distance_int(struct city* city1, struct city* city2){

    int x_distance = abs(city1->x_coordinate - city2->x_coordinate);
    int y_distance = abs(city1->y_coordinate - city2->y_coordinate);
    double distance = sqrt( (x_distance*x_distance) + (y_distance*y_distance) );

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

double** create_distance_matrix_int(struct city* cities, int dimension){
    double** distance_m = malloc(dimension * sizeof(double*));

    for (int i = 0; i < dimension; i++) {
        distance_m[i] = malloc(dimension * sizeof(double));

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

void initialize_with_nearest_neighbor(struct population* pop, struct city* cities,float** distance_m,int dimension, double percentage){

    int create_n_random = pop->population_count * percentage;
    //srand(time(NULL));
    int random = 0;
    for(int i = 25; i< 25 + create_n_random; ++i){
        random = get_one_random_number(dimension);
        //Put first random city in the list
        //pop->individuals[i].cities[0] = cities[random];
        memcpy(&pop->individuals[i].cities[0], &cities[random], sizeof(struct city));

        node_t * head = NULL;
        head = malloc(sizeof(node_t));
        head->val = random;
        head->next = NULL;
        node_t* last = head;
        for(int j = 1; j< 280 ;j++){
            int nearest_neighbor_distance= BIGINT;
            int nearest_neighbor=0;
            for(int x = 0; x < dimension; ++x){
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
            //pop->individuals[i].cities[j] = cities[nearest_neighbor];
            memcpy(&pop->individuals[i].cities[j], &cities[nearest_neighbor], sizeof(struct city));

            random = nearest_neighbor;
            node_t* next = malloc(sizeof(node_t));
            next->val = nearest_neighbor;
            last->next = next;
            next->next = NULL;
            last = next;

        }


        node_t * tmp;

        while (head != NULL) {
            tmp = head;
            head = head->next;
            free(tmp);

        }
    }

}

void initialize_randomly(struct population* pop, struct city* cities, double percentage){

    int create_n_random = pop->population_count * percentage;
    //srand(time(NULL));

    for(int i = 0; i< create_n_random; i++){
        //int random = rand() % 280;
        int* r = malloc(280 * sizeof(int));

        for(int i=0;i<280;++i){
            r[i]=i;
        }

        for (int i = 280-1; i >= 0; --i){
            //generate a random number [0, n-1]
            int j = get_one_random_number(i + 1);

            //swap the last element with element at random index
            int temp = r[i];
            r[i] = r[j];
            r[j] = temp;
        }

        for(int j = 0; j< 280 ;j++){
            memcpy(&pop->individuals[i].cities[j], &cities[r[j]], sizeof(struct city));
            //pop->individuals[i].cities[j] = cities[r[j]];//TODO check this seems not correct

        }
        free(r);
    }
}

float calculate_fitness(struct city* cities, float** distance_m, int dimension){
    float fitness_value = 0;
    for(int i = 0;i<dimension - 1;i++){
        fitness_value = fitness_value + distance_m[cities[i].city_id - 1][cities[i+1].city_id - 1];
       // printf("fitness value: %f", fitness_value);
    }
    fitness_value = fitness_value + distance_m[cities[dimension - 1].city_id - 1][cities[0].city_id - 1];

    return fitness_value;
}
double calculate_fitness_int(struct city* cities, double** distance_m, int dimension){
    double fitness_value = 0;
    for(int i = 0;i<dimension - 1;i++){
        fitness_value = fitness_value + distance_m[cities[i].city_id - 1][cities[i+1].city_id - 1];
        // printf("fitness value: %f", fitness_value);
    }
    fitness_value = fitness_value + distance_m[cities[dimension - 1].city_id - 1][cities[0].city_id - 1];

    return fitness_value;
}

struct city** perform_OX(struct city *parent1, struct city *parent2, int dimension) {

    struct city* child1 = calloc(dimension, sizeof(struct city));
    struct city* child2 = calloc(dimension, sizeof(struct city));

    int crossover_start = get_one_random_number(dimension);

    int crossover_end = get_one_random_number(dimension - crossover_start) + crossover_start;


//    memset(child1,0,dimension* sizeof(struct city));
//    memset(child2,0,dimension* sizeof(struct city));

    //First copy subset from parents to children
    for(int i = crossover_start;i<crossover_end;++i){
        memcpy(&child1[i], &parent2[i], sizeof(struct city));
        memcpy(&child2[i], &parent1[i], sizeof(struct city));
//        child1[i] = parent2[i];
//        child2[i] = parent1[i];
    }


    for(int i = crossover_end; i< dimension; ++i){
        int exists_in_subset = 0;
        for(int i1 = crossover_start; i1 < crossover_end; ++i1){
            if(parent2[i].city_id == child2[i1].city_id)
                exists_in_subset = 1;
        }
        if(0 == exists_in_subset)
            memcpy(&child2[i], &parent2[i], sizeof(struct city));

        exists_in_subset = 0;
        for(int i1 = crossover_start; i1 < crossover_end; ++i1){
            if(parent1[i].city_id == child1[i1].city_id)
                exists_in_subset = 1;
        }
        if(0 == exists_in_subset)
            memcpy(&child1[i], &parent1[i], sizeof(struct city));

    }

    int parent1_index = 0;
    int parent2_index = 0;

    for(int i = 0; i<dimension;++i){
        if(0 == child2[i].city_id){
            for(int j = parent2_index;j<dimension;++j){
                int exists_in_subset = 0;
                for(int i1 = 0; i1 < dimension; ++i1){
                    if(parent2[j].city_id == child2[i1].city_id)
                        exists_in_subset = 1;
                }
                if(0 == exists_in_subset){
                    memcpy(&child2[i], &parent2[j], sizeof(struct city));
                    //child2[i] = parent2[j];
                    parent2_index = j;
                    break;
                }
            }
        }

    }

    for(int i = 0; i<dimension;++i){
        if(0 == child1[i].city_id){
            for(int j = parent1_index;j<dimension;++j){
                int exists_in_subset = 0;
                for(int i1 = 0; i1 < dimension; ++i1){
                    if(parent1[j].city_id == child1[i1].city_id)
                        exists_in_subset = 1;
                }
                if(0 == exists_in_subset){
                    memcpy(&child1[i], &parent1[j], sizeof(struct city));
                    //child1[i] = parent1[j];
                    parent1_index = j;
                    break;
                }
            }
        }

    }

//    struct offsprings* offsprings = malloc(sizeof(struct offsprings));
//    offsprings->offspring1 = malloc(280* sizeof(struct city));
//    offsprings->offspring2 = malloc(280* sizeof(struct city));
//
//    offsprings->offspring1 = child1;
//    offsprings->offspring2 = child2;
     struct city** childeren = malloc(2 * sizeof(struct city*));
    childeren[0] = calloc(dimension, sizeof(struct city));
    childeren[1] = calloc(dimension, sizeof(struct city));
//    memset(childeren[0],0,dimension* sizeof(struct city));
//    memset(childeren[1],0,dimension* sizeof(struct city));
    memcpy(&childeren[0][0], child1,  dimension *sizeof(struct city));
    memcpy(&childeren[1][0], child2,  dimension *sizeof(struct city));
    free(child1);
    free(child2);
//
    return childeren;
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


void perform_tournament(struct population* pop,size_t tournament_size, size_t dimension, float** distance_m, int* picked_parent_index){

    int random_individual = 0;
    struct individual tournament[tournament_size];
    int picked_numbers[tournament_size];
    int* random_numbers = malloc(tournament_size * sizeof(int));
    get_n_different_random_number(random_numbers, tournament_size, pop->population_count);

    float min_fitness = 99999;
    int current_individual = 0;
    for (int i = 0; i < tournament_size; ++i) {
        float individual_fitness = calculate_fitness(pop->individuals[random_numbers[i]].cities, distance_m, dimension);
        if(individual_fitness < min_fitness){
            min_fitness = individual_fitness;
            current_individual = random_numbers[i];
        }

    }
    free(random_numbers);
    *picked_parent_index = current_individual;
   // return pop->individuals[current_individual].cities;
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
struct city** perform_SCX(struct city* parent1, struct city* parent2, float** distance_m, int dimension){

    int size = 0;
    int index = 0;
    int cost_to_pick_from_p2 = 0;
    int cost_to_pick_from_p1 = 0;
    struct city* offspring = calloc(dimension , sizeof(struct city));
    //memset(offspring,0,dimension* sizeof(struct city));
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

    struct city** childeren = malloc(1 * sizeof(struct city));
    childeren[0] = calloc(dimension, sizeof(struct city));
    //memset(childeren[0],0,dimension* sizeof(struct city));
    memcpy(&childeren[0][0], offspring,  dimension *sizeof(struct city));
    free(offspring);

    return childeren;

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

    /*first_random_gene = 1;
    second_random_gene = 4;*/
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

    if(*random1 < *random2){
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

int find_worst_individual_index(struct population* pop, float** distance_m, size_t dimension){
    float worst = 0;
    int index = 0;
    for (int i = 0; i < 50; ++i) {
        int fit = calculate_fitness(pop->individuals[i].cities, distance_m, dimension);
        if(fit > worst){
            worst = fit;
            index = i;
        }


    }

    return index;
}

void replace_with_worst_individual(struct population* pop, struct city* child, float** distance_m, size_t dimension){
    int worst = find_worst_individual_index(pop, distance_m, dimension);

    memcpy(pop->individuals[worst].cities, child,  dimension *sizeof(struct city));

}

void replace_with_random_individual(struct population* pop, struct city* child, size_t dimension){
    int random = get_one_random_number(pop->population_count);

    memcpy(pop->individuals[random].cities, child,  dimension *sizeof(struct city));
}

int find_best_individual(struct population* pop,float** distance_m, size_t dimension){

    float best = 99999999;
    int index = 0;
    for (int i = 0; i < pop->population_count; ++i) {
        int fit = calculate_fitness(pop->individuals[i].cities, distance_m, dimension);
        if(fit < best){
            best = fit;
            index = i;
        }


    }

    return index;
}

void two_opt(struct city* individual,float** distance_m, int dimension){
    int first_index = get_one_random_number(dimension);
    int second_index = get_one_random_number(dimension - first_index);
    second_index = first_index + second_index;
//    if(*random_1 < random_2){
//        first_index = *random_1;
//        second_index = *random_2;
//    }
//    else{
//        first_index = *random_2;
//        second_index = *random_1;
//    }
    first_index = 2;
    second_index = 7;

    struct city* copy_of_individual = calloc(dimension, sizeof(struct city));
    //memset(copy_of_individual, 0 , dimension* sizeof(struct city));
    memcpy(copy_of_individual, individual, dimension * sizeof(struct city));


    int intermediate_1_size = second_index - first_index - 1;

    struct city* intermediate = malloc(intermediate_1_size* sizeof(struct city));

    for (int i = 0; i < intermediate_1_size; ++i) {
        memcpy(intermediate + i, copy_of_individual + second_index - i - 1, (sizeof(struct city)));

    }

    int intermediate_2_size = dimension - second_index - 1;

    struct city* intermediate2 = malloc(intermediate_2_size* sizeof(struct city));

    memcpy(intermediate2, copy_of_individual + second_index + 1, (intermediate_2_size * sizeof(struct city)));

    memmove(copy_of_individual + first_index + 1, copy_of_individual + second_index, sizeof(struct city));
    memmove(copy_of_individual + first_index + 2, intermediate, intermediate_1_size * sizeof(struct city));
    memmove(copy_of_individual + first_index + 2 + intermediate_1_size, intermediate2, intermediate_2_size * sizeof(struct city));

    int fitness_of_individual = calculate_fitness(individual, distance_m, dimension);
    int fitness_of_new_tour = calculate_fitness(copy_of_individual, distance_m, dimension);

    if(fitness_of_new_tour < fitness_of_individual)
        memcpy(individual, copy_of_individual, dimension * sizeof(struct city));

}

void perform_2opt(struct population* pop, float** distance_m, size_t dimention, size_t n_times, size_t m_random_ind){

    int best_individual = find_best_individual(pop, distance_m, dimention);
    for (int k = 0; k < n_times; ++k) {
        two_opt(pop->individuals[best_individual].cities,distance_m, dimention);

    }

    int random_individual = 0;
    for (int i = 0; i < m_random_ind; ++i) {
        do{
            random_individual = get_one_random_number(dimention);
        }while(random_individual == best_individual);

        for (int j = 0; j < n_times; ++j) {
            two_opt(pop->individuals[random_individual].cities,distance_m, dimention);
        }
    }

}

float calculate_average_fitness(struct population* pop, float** distance_m, size_t dimension){
    float fitness = 0;
    for (int i = 0; i < pop->population_count; ++i) {
        fitness = fitness + calculate_fitness(pop->individuals[i].cities, distance_m, dimension);
    }

    return fitness / pop->population_count;


}

int main(int argc, char *argv[]){
    printf("------------------------------\n");
    printf("Experiment start\n");
    printf("------------------------------\n");


    //Parameters are set via arg
    int iteration_count = 0, tournament_size = 0,  show_always = 0;
    enum CrossOverType xover_type;
    enum MutationType mutation_type;

    //variables to read from file
    FILE *fp;
    char line[200];
    struct city* cities;
    fp = fopen("../a280.tsp", "r"); // read mode

    if (fp == NULL)
    {
        perror("Error while opening the file.\n");
        exit(EXIT_FAILURE);
    }

    //Read arguments for algorithm
    if(argc == 6){
        if(argv[1] != NULL){
            sscanf(argv[1], "%d", &iteration_count);
        }
        if(argv[2] != NULL){
            sscanf(argv[2], "%d", &xover_type);
        }
        if(argv[3] != NULL){
            sscanf(argv[3], "%d", &mutation_type);
        }
        if(argv[4] != NULL){
            sscanf(argv[4], "%d", &tournament_size);
        }
        if(argv[5] != NULL){
            sscanf(argv[5], "%d", &show_always);
        }
    } else
        return 0;


    printf("/-----------------------------------/\n");
    printf("Arguments read from command line\n");
    printf("/-----------------------------------/\n");


    //Distance Matrix
    float** distance_m;

    //Random mutation
    if(Random == mutation_type){
        int rand = get_one_random_number(3);
        mutation_type = rand;
    }

    const int mutation_rate = 0;
//    //Crossover function pointer
//    struct city* (*cross_over_function_ptr)(int);
//    switch(xover_type){
//        case  OX:
//
//    }

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
                     cities = calloc( dimension , sizeof(struct city));

                }
            }
            continue;
        }

        struct city *n= calloc(1,sizeof(struct city));

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
    double** dist_int = create_distance_matrix_int(cities, dimension);
    //printf("dimension between 1 and 2 is: %f", distance_m[0][1]);
    //Initialize a population;
    struct population* pop = calloc(1, sizeof(struct population));
    pop->population_count = 50;
    pop->individuals = calloc(50 , sizeof(struct individual));
    for(int i = 0; i < pop->population_count; i++){
        pop->individuals[i].cities = calloc(dimension , sizeof(struct city));
    }

    //Initialize half of the population randomly
    initialize_randomly(pop, cities, 0.5);

    //struct city* tour = (struct city*)malloc( 280 * sizeof(struct city));
    initialize_with_nearest_neighbor(pop, cities,distance_m, dimension, 0.5);

    struct city* child1;
    struct city* child2;
    int child_count = 0;
    if(OX == xover_type){

         child_count = 2;
    }else if(SCX == xover_type){

         child_count = 1;
    }
    struct city* parent1;
    struct city* parent2;

//    struct city* cities2 = malloc(7 * sizeof(struct city));
//    cities2[0].city_id = 1;
//    cities2[1].city_id = 6;
//    cities2[2].city_id = 2;
//    cities2[3].city_id = 4;
//    cities2[4].city_id = 3;
//    cities2[5].city_id = 5;
//    cities2[6].city_id = 7;
//    struct city** off1 = perform_OX(cities1, cities2, 7);
    fp = fopen("../a280.opt.tour", "r"); // read mode
    struct city* opt = calloc( 280 , sizeof(struct city));

    if (fp == NULL)
    {
        perror("Error while opening the file.\n");
        exit(EXIT_FAILURE);
    }
    pass_first_lines = 0;
    i = 0;
    while ( ( fgets ( line, sizeof ( line), fp)) != NULL) {
        if(pass_first_lines++ < 4){
            continue;
        }

        struct city *n= calloc(1,sizeof(struct city));

        if ( ( sscanf ( line, "%d"
                , &n->city_id )) == 1){
            n->y_coordinate = 0;
            n->x_coordinate = 0;
            opt[i] = *n;
            i++;
        }
    }
    fclose(fp);
    double asf = calculate_fitness_int(opt, dist_int, 280);
    FILE *f = fopen("../a280tsp.txt", "a");

    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    int best_before_experiment = find_best_individual(pop, distance_m, dimension);
    int best_tour_before_exp = calculate_fitness(pop->individuals[best_tour_before_exp].cities, distance_m, dimension);
    fprintf(f, "Best fitness before exp: %d \n", best_tour_before_exp);
    fprintf(f, "Experiment with XOver: %d, mutation: %d, iteration: %d \n", xover_type, mutation_type, iteration_count);

    int best_result = 0;
    int average_result = 0;
    int iter = 1;

    while(iter <= iteration_count){
        int* index_parent1 = malloc(sizeof(int));
        int* index_parent2 = malloc(sizeof(int));

        do{
            perform_tournament(pop, tournament_size, dimension, distance_m, index_parent1);
            perform_tournament(pop, tournament_size, dimension, distance_m, index_parent2);
        }while(*index_parent1 == *index_parent2);

       struct city** off;
       if(OX == xover_type){

           off = perform_OX(pop->individuals[*index_parent1].cities, pop->individuals[*index_parent2].cities, dimension);
       }else if(SCX == xover_type){
           off = perform_SCX(pop->individuals[*index_parent1].cities, pop->individuals[*index_parent2].cities, distance_m, dimension);
       }

       free(index_parent1);
       free(index_parent2);
       int mutation = get_one_random_number(9);
       if(mutation <= mutation_rate){
           printf("/-----------------------------/\n");
           printf("mutation  iter: %d\n", iter);
           printf("/-----------------------------/\n");
           if(ISM == mutation_type){
               for (int j = 0; j < child_count; ++j) {
                   perform_Insert_mutation(&off[j][0], dimension);
               }
           }else if(IVM == mutation_type){
               for (int j = 0; j < child_count; ++j) {
                   perform_inversion_mutation(&off[j][0], dimension);
               }
           }else if(SM == mutation_type){
               for (int j = 0; j < child_count; ++j) {
                   perform_swap_mutation(&off[j][0], dimension);
               }
           }
       }
       int random_child = get_one_random_number(child_count);
       replace_with_worst_individual(pop, off[random_child], distance_m, dimension);
       if(child_count > 1)
           replace_with_random_individual(pop, off[1 - random_child], dimension);

       free(off[0]);
       if(child_count > 1)
        free(off[1]);

       if(1 == show_always || (iter == 1000 || iter == 5000 || iter == 10000)){
            int current_best = find_best_individual(pop, distance_m, dimension);
            float current_best_fitness = calculate_fitness(pop->individuals[current_best].cities, distance_m, dimension);
            float average = calculate_average_fitness(pop, distance_m, dimension);
            fprintf(f, "/---------------------------------- /\n");
            fprintf(f, "Iteration: %d \n", iter);
            fprintf(f, "Best Tour: %f \n", current_best_fitness);
            fprintf(f, "Average fitness: %f \n", average);
            fprintf(f, "/---------------------------------- /\n");
       }

       iter++;
    }

    int best = find_best_individual(pop, distance_m, dimension);
    float best_tour = calculate_fitness(pop->individuals[best].cities, distance_m, dimension);

//    for (int i1 = 0; i1 < dimension; ++i1) {
//        fprintf(f, "%d\n", pop->individuals[best].cities[i1].city_id);
//
//    }
    fprintf(f, "/---------------------------------- /\n");
    fprintf(f, "Experiment Finished\n");
    fprintf(f, "Best tour: %f \n", best_tour);
    fprintf(f, "/---------------------------------- /\n");
    fclose(f);

//    struct city* cities3 = malloc(9 * sizeof(struct city));
//    cities3[0].city_id = 1;
//    cities3[1].city_id = 2;
//    cities3[2].city_id = 3;
//    cities3[3].city_id = 4;
//    cities3[4].city_id = 5;
//    cities3[5].city_id = 6;
//    cities3[6].city_id = 7;
//    cities3[7].city_id = 8;
//    cities3[8].city_id = 9;


//    perform_inversion_mutation(cities3, 9);
//

        //int random = get_one(7);
        //int* random = get_one_random_number(7);

    free(distance_m);
    for(int i = 0; i < pop->population_count; i++){
        free(pop->individuals[i].cities);
    }
    free(pop->individuals);
    free(pop);
    return 0;


}
