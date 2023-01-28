#include "network.h"

void print_usage0(){
    printf("===========================\n");
    printf("Author:   ChiLianjiang\n");
    printf("E-mail:   chilianjiang@126.com\n");
    printf("Date:     2021-03-24\n");
    printf("Version:  1.0\n\n");
    printf("fastHapNetwork original_tcs [arguments]\n");
    printf("fastHapNetwork modified_tcs [arguments]\n");
    printf("fastHapNetwork msn          [arguments]\n");
    printf("fastHapNetwork mjn          [arguments]\n\n");
    exit(0);
}

int main(int argc, char **argv){
    if(argc<2) print_usage0();
    argc--;
    argv++;

    if(strcmp(argv[0], "original_tcs")==0) main_tcs(argc, argv);
    else if(strcmp(argv[0], "modified_tcs")==0 || strcmp(argv[0], "msn")==0 || strcmp(argv[0], "mjn")==0) main_mjn(argc, argv);
    else print_usage0();

    return 0;
}







