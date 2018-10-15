#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

const int BIGINT = 9999999;
enum CrossOverType {OX = 0, SCX = 1};
enum MutationType {ISM = 0, IVM = 1, SM = 2, Random = 3};

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



//void initialize_with_nearest_neighbor(struct population* pop, struct city* cities, double percentage){
//
//    int create_n_random = pop->population_count * percentage;
//    srand(time(NULL));
//
//    for(int i = 0; i< create_n_random; i++){
//        int random = rand() % 280;
//
//        for(int j = 0; j< 270 ;j++){
//            pop->individuals[i].cities[j] = cities[r[j]];
//
//        }
//
//    }
//
//
//}

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
    while ( ( fgets ( line, sizeof ( line), fp))) {
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
    printf("dimension between 1 and 2 is: %f", distance_m[0][1]);
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
    //initialize_with_nearest_neighbor(pop, cities, 0.5);




    return 0;


    printf("Hello, World!\n");
    return 0;


}