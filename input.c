#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "input.h"
#include "xmmintrin.h"


#define INF 9999999

static int *in_matrix;// __attribute__((aligned (32)));

static void find32(int val){

    while(val % 32 != 0){
        val++;
    }
    printf("/32 is %d\n",val);
}

static void usage(const char *pname){

    printf("Usage: %s [OPTION]...\n"
           "\n"
           "  -n NUM     Set total number of nodes in the graph.\n"
           "  -e NUM     Set total number of edges in the graph.\n"
           "  -c FILE    Read input from FILE.\n"
           "  -p NUM     Number of threads to use (when applicable).\n"
           "  -r NUM     Report general statistic at the end.\n"
           "  -m NUM     Print the resulting matrix.\n"
           "  -v NUM     Validates results against sequential version (slows down execution)\n"
           "  -h         Print this help.\n"
           ,pname);
    exit(0);
}


/*This function is specific for V#-E## files*/
static void read_input_file(const char *fname, size_t N, size_t E, int *data, struct parameters* p){
    FILE *f;
    int i,j;
    int current_line=0;
    int line_index=0;
    int weight = 0;
    int max_weight = 0;
    int test=100;
    char buff[256]; //Seems to be ok
    int temp_n, temp_e;
    /*Init data structure to INT_MAX. */
    for(int m=0; m<N; ++m){
        for(int k=0; k<N; ++k){
            data[m*N+k] = (m==k)? 0 : INF;
        }
    }


    printf("Reading input file %s\n", fname);
    if (!(f = fopen(fname, "r"))){printf("Could not open input file, exiting...\n");exit(1);}

    while(fgets(buff,256,f)!=NULL){
        int field;
        int n;
        char *begin = buff;
        /*Read first 2 lines:(Nodes and Edges)*/
        current_line++;
        if(current_line==1) temp_n = atoi(buff);
        if(current_line==2) temp_e = atoi(buff);

        /*Use sscanf to split lines into values*/
        if(current_line>2){
            if(temp_n!=p->n || temp_e!=p->e){
                printf("Either N or E provided through args are not equal to the values specified in input file. Exiting...\n");
                exit(1);
            }
            while (sscanf(begin, "%d%n", &field, &n)==1){
                line_index++;
                if(line_index==1) i = field-1; //from
                if(line_index==2) j = field-1; //to
                if(line_index==3) weight = field; //weight

                begin+=n; //add offset to starting position
            }
            //After line is read, store the input.
            data[i * N + j] = (i==j)? 0 : weight; //Negates case where input is 2 2 515, meaning path to self
            {i=0; j=0; weight=0;}       //reset
        }
        line_index=0;
    }

    fclose(f);
    if(temp_n!=p->n || temp_e!=p->e){
        printf("Either N or E provided through args are not equal to the values specified in input file. Exiting...\n");
        exit(1);
    }
}



void read_parameters(struct parameters* p, int argc, char **argv){
    const char * input_fname=0;
    int ch;

    /*Defaults- Keep the debug file in the same dir*/
    p->n = 5;
    p->e = 10;
    p->nthreads = 1;
    p->print_results = 1;
    p->print_output = 0;
    p->validate = 0;
    input_fname = "V5-E10";


    while ((ch = getopt(argc, argv, "c:e:E:hH:n:N:p:m:r:v")) != -1)
    {
        switch(ch) {
        case 'c': input_fname = optarg; break;
        case 'n': case 'N': p->n = strtol(optarg, 0, 10); break;
        case 'e': case 'E': p->e = strtol(optarg, 0, 10); break;
        case 'p': p->nthreads = strtol(optarg, 0, 10); break;
        case 'r': p->print_results = 1; break;
        case 'm': p->print_output = 1; break;
        case 'v': p->validate = 1; break;
        case 'h': default: usage(argv[0]);
        }
    }


    printf("Parameters:\n"
        "  -n %zu # number of nodes\n"
        "  -e %zu # number of edges\n"
        "  -c %s # input file for graph data\n"
        "  -p %zu # number of threads (if applicable)\n"
        "  -r %zu # print execution results\n"
        "  -m %zu # print resulting matrix\n"
        "  -v %d #validates parallel version\n",
        p->n, p->e, 
        input_fname ? input_fname : "(none)",
        p->nthreads, p->print_results, p->print_output, p->validate);

    if(!p->n || !p->e){printf("n or e werer empty during init, exiting...\n"); exit(1);}

    int size = p->n * p->n;
    //Round up to nearest multiple of 32 if not already a multiple of 32. Useful for allignment
    int elements = size-1;
    elements = elements >>5;
    elements+=1;
    elements = elements <<5;
    /*Use _mm_malloc for the SIMD versions*/
    in_matrix = _mm_malloc(elements * sizeof(int),32); 
    /*Use this for seq*/
    //in_matrix = malloc(size * sizeof(int));

    if(in_matrix==NULL){
        printf("Malloc failed for input matrix, exiting...\n");
        exit(1);
    } else {
        read_input_file(input_fname, p->n, p->e, in_matrix, p);
        p->input_matrix = in_matrix;
    }
}