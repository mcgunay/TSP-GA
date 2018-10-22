//
// Created by mcangny on 13/10/18.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(){

    printf("experiments start");
    int status;
//    for (int i = 0; i < 100; ++i) {
//        status = system("./test1 10000 0 0 5 0");
//    }
//    for (int i = 0; i < 100; ++i) {
//        status = system("./test1 10000 0 1 5 0");
//    }
//    for (int i = 0; i < 100; ++i) {
//        status = system("./test1 10000 0 2 5 0");
//    }
//    for (int i = 0; i < 100; ++i) {
//        status = system("./test1 10000 0 3 5 0");
//    }
//    for (int i = 0; i < 100; ++i) {
//        status = system("./test1 10000 1 0 5 0");
//    }
//    for (int i = 0; i < 100; ++i) {
//        status = system("./test1 10000 1 1 5 0");
//    }
//    for (int i = 0; i < 100; ++i) {
//        status = system("./test1 10000 1 2 5 0");
//    }
//    for (int i = 0; i < 100; ++i) {
//        status = system("./test1 10000 1 3 5 0");
//    }
     status = system("./test1 10000 0 0 5 0");
     status = system("./test1 10000 0 1 5 0");
     status = system("./test1 10000 0 2 5 0");
     status = system("./test1 10000 0 3 5 0");
     status = system("./test1 10000 1 0 5 0");
     status = system("./test1 10000 1 1 5 0");
     status = system("./test1 10000 1 2 5 0");
     status = system("./test1 10000 1 3 5 0");


    printf("status: %d", status);
    printf("experiments end");

    return 0;
}