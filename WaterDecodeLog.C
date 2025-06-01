#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

#define max(a,b) \
({ __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })


double *F_store;
double *B_store;
double drift;


// Define a struct to hold necessary data
typedef struct {
    double Pi, Pd, Pt;
    double Ps;
    int N, Nr, Nd;
    int *info;
    int *r, *w, *d;
    double *priors;
    int *data; // only to display success rate
} ChannelData;

double forward(ChannelData *data, int i, int j);
double backward(ChannelData *data, int i, int j);


bool over_limit(int i, int j){
    int a  = abs(i - j);
    if (a < 10){
        return false;
    }
    float max_drift = 5 * pow(drift *i, 0.5);
    if (a > max_drift){
        return true;
    }
    else{
        return false;
    }
}

void initialise_channel_data(ChannelData *data) {
    int s = (data->N + 1) * (data->Nr + 1);
    for (int i = 0; i < s; i++) {
        F_store[i] = 1.0;
        B_store[i] = 1.0;
    }
}

int *int_read_from_file(const char *filename, int size) {
    // fprintf(stderr, "Reading from file %s size %d \n", filename, size);
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file: %s\n", filename);
        exit(1);
    }
    int *array = (int *)malloc(size * sizeof(int));
    int c;
    for (int i = 0; i < size; i++) {
        c = fscanf(file, "%d", &array[i]);
    if (c==EOF) 
      { if (i>0)
        { fprintf(stderr,
          "Warning: Short block (%d long) in %s \n",i, &filename);
        }
      }
    }

    fclose(file);
    return array;
}

double logaddexp(double x1, double x2){
    if (x1 == -INFINITY and x2 == -INFINITY){
        return -INFINITY;
    }
    double a = max(x1, x2);
    double b;
    if (a > x1){
        b = x1;
    }
    else{
        b = x2;
    }
    double l = log(1 + exp(- (a - b)));
    // fprintf(stderr, "Log adding %f to %f  and got (%f)", x1, x2, a + l);
    return a + l;
}

double *d_read_from_file(const char *filename, int size) {
    // fprintf(stderr, "Reading from file %s with size %d \n", filename, size);
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file: %s\n", filename);
        exit(1);
    }
    int c;
    double *array = (double *)malloc(size * sizeof(double));

    for (int i = 0; i < size; i++) {
        fscanf(file, "%lf", &array[i]);
    }
    return array;

    fclose(file);
}


double gamm(ChannelData *data, int i1, int j1, int i2, int j2) {
    double g;
    if (i2 == i1 && j2 == j1 + 1) {
        g = data->Pi / 2;
    } else if (i2 == i1 + 1 && j2 == j1) {
        g = data->Pd;
    } else if (i2 == i1 + 1 && j2 == j1 + 1) {
        bool found = false; // if bit is data bit
        int ind;
            for (int k=0; k<data->Nd; k++) 
            { 
            if (data->info[k] == i1) 
            { 
                found = true;
                ind = k;
                break; 
            } 
            } 
        if (!found) {
            g = (data->r[j1] == data->w[i1]) ? data->Pt * (1 - data->Ps) : data->Pt * data->Ps;
        } else {
            double p = (1.0 - data->priors[ind]); //probability that bit is 0
            g = ((data->r[j1] + data->w[i1]) % 2 == 1) ? 
                (p * data->Pt * data->Ps + (1 - p) * data->Pt * (1 - data->Ps)) :
                (p * data->Pt * (1 - data->Ps) + (1 - p) * data->Pt * data->Ps);
        }
    }
    if (g > 1){
        fprintf(stderr, "Got large gamma value of %f", g);
    }
    double lg;
    
    if(g == 0){
        lg =  -INFINITY;
    }
    else{
        lg =  log(g);
    }
    // fprintf(stderr, "Got gamma value of %12.5e for %d %d -> %d %d \n", log(g), i1, j1, i2, j2);
    return lg;
}

double compute_forward(ChannelData *data, int i, int j) {
    // fprintf(stderr, "Computing forwards for %d %d \n", i, j);
    double f;
    if (i > data->N || j > data->Nr || i < 0 || j < 0) {
        fprintf(stderr, "Outside\n" );
        f = -INFINITY;
    } else if (i == 0 && j == 0) {
        // fprintf(stderr, "Boundary\n" );
        f  =  0.0;
    } else if (i == 0) {
        f = forward(data, i, j - 1) + gamm(data, i, j - 1, i, j);
    } else if (j == 0) {
        f = forward(data, i - 1, j) + gamm(data, i - 1, j, i, j);
    } else {
        double f1 = forward(data, i - 1, j) + gamm(data, i - 1, j, i, j);
        double f2 = forward(data, i, j - 1) + gamm(data, i, j - 1, i, j);
        double f3 = forward(data, i - 1, j - 1) + gamm(data, i - 1, j - 1, i, j);
        f = logaddexp(f1, logaddexp(f2, f3));
    }
    // fprintf(stderr, "Got forward value of %12.5e for %d %d \n", f, i, j);

    F_store[data->N*j + i] = f;
    return f;
}

double compute_backward(ChannelData *data, int i, int j) {
    //fprintf(stderr, "Computing backwards for %d %d \n", i, j);
    double b;
    if (i == data->N && j == data->Nr) {
        // fprintf(stderr, "Hit Boundary at %d %d\n", i, j );
        b = 0;
    } else if (i == data->N) {
        b = -INFINITY;
    } else if (j == data->Nr) {
        b = backward(data, i + 1, j) + gamm(data, i, j, i + 1, j);
    } else if (i > data->N || j > data->Nr || i < 0 || j < 0) {
        b = -INFINITY;
    } else {
        double b1 = backward(data, i + 1, j) + gamm(data, i, j, i + 1, j);
        // fprintf(stderr, "Got b1 as %f", b1);
        double b2 = backward(data, i, j + 1) + gamm(data, i, j, i, j + 1);
        // fprintf(stderr, "Got b2 as %f", b1);
        double b3 = backward(data, i + 1, j + 1) + gamm(data, i, j, i + 1, j + 1);
        // fprintf(stderr, "Got b3 as %f", b1);
        b = logaddexp(b1, logaddexp(b2, b3));
    }
    B_store[data->N*j + i] = b;
    if (i >= data->N - 2 || j >= data->Nr - 2) {
        // fprintf(stderr, "Got backward value of %12.5e for %d %d \n", b, i, j);
    }
    // fprintf(stderr, "Got backward value of %12.5e for %d %d \n", b, i, j);
    if (b != b){
        exit(1);
    }
    return b;
}

double forward(ChannelData *data, int i, int j) {
    if (F_store[data->N*j + i] != 1.0) {
        // printf("Got f value from store");
        return F_store[data->N*j + i];
    }
    if (over_limit(i, j)){
        F_store[data->N*j + i] = -INFINITY;
        return -INFINITY;
    }
    return compute_forward(data, i, j);
}

double backward(ChannelData *data, int i, int j) {
    if (B_store[(data->N * j) + i] != 1.0) {
        // fprintf(stderr, "Got from store for %d %d value is %f \n", i, j, B_store[(data->N * j) + i]);
        return B_store[(data->N * j) + i];
    }
    if (over_limit(i, j)){
        B_store[data->N*j + i] = -INFINITY;
        // fprintf(stderr, "OL");
        return -INFINITY;
    }
    // fprintf(stderr, "Computing");
    return compute_backward(data, i, j);
}

double likelihood(ChannelData *data, int i, int val){
    // fprintf(stderr, "Attempting likelihood now \n");
    double l = -INFINITY;
    for (int j = 0; j < data->Nr; j++){
        if (over_limit(i, j)){
            continue;
        }
        double diag;
        if (data->r[j] == (val + data->w[i]) % 2){
            diag = log(data->Pt * (1 - data->Ps));
        }
        else{
            diag = -INFINITY; // as ps = 0 at the moment
            //diag = data->Pt*data->Ps;
        }
        double l2 = forward(data, i, j) + diag + backward(data, i + 1, j + 1);
        l = logaddexp(l, l2);

        double l3 = forward(data, i, j) + log(data->Pd) + backward(data, i + 1, j);
        l = logaddexp(l, l3);
    }
    double l4 = forward(data, i, data->Nr) + log(data->Pd) + backward(data, i + 1, data->Nr);
    l = logaddexp(l, l4);

    return l;
    
}

int decode(ChannelData *data, FILE* llr_file){
    // fprintf(stderr, "Attempting decode now");
    double *llrs = (double*)malloc(data->Nd * sizeof(double));
     if (llrs==0)
    { fprintf(stderr,"Ran out of memory (while trying to allocate bytes)\n");
        exit(1);
    }
    int successes = 0;
    for (int j = 0; j < data->Nd; j++ ){
        double pr[2];
        double p_win = -INFINITY;
        int ans;
        for (int v = 0; v < 2; v ++){
            pr[v] = likelihood(data, data->info[j], v);
            if (pr[v] > p_win){
                p_win = pr[v];
                ans = v;
            }
            // fprintf(stderr, "Likelihood: %.10f ", pr[v]);
            // if (pr[v] == 0){
            //     ;
            //    printf("Got zero p! for data bit %d value %d \n", j, v);
            // }
        }
        if (data->d[j] == ans){
            successes = successes + 1;
        }
    
        llrs[j] = pr[0] -  pr[1];
        //printf("LLR %d: %f  \n", j, llrs[j]);
        fprintf(llr_file, "%.6f ", llrs[j]);
            
        
    }
    printf("%d/%d successes \n", successes, data->Nd);
    free(F_store);
    free(B_store);
    return 0;
}

int main(int argc, char *argv[]) {
   if (argc != 11) {
        fprintf(stderr, "Usage: %s <priors_file> <llrs_file> <received_file> <info_file> <w_file> <N> <Nr> <Nd> <Pipd> <d_file>\n", argv[0]);
        return 1;
    }
    ChannelData data;
    data.N = atoi(argv[6]);
    // fprintf(stderr, "Got N value %d \n", data.N);
    data.Nr = atoi(argv[7]);
    // fprintf(stderr, "Got Nr value %d \n", data.Nr);

    F_store = (double*)malloc((data.N + 1) * (data.Nr + 1) * sizeof(double));
    B_store = (double*)malloc((data.N + 1) * (data.Nr + 1) * sizeof(double));

    initialise_channel_data(&data); 
    char *llr_file = argv[2];
    // printf("Writing llrs to %s", llr_file);
    FILE *llrf;

    llrf = fopen(llr_file, "w");
    if (llrf==NULL)
    { fprintf(stderr,"Can't create file for decoded data: %s\n",llr_file);
        exit(1);
    }
    
    // setup currently that pi = pd
    data.Nd = atoi(argv[8]);
    data.Pi = atof(argv[9]);
    data.Pd = atof(argv[9]);
    drift = (data.Pi + data.Pd) / (1 - data.Pi); //if pi = pd only
    data.Ps = 0;
    data.Pt = 1.0 - (data.Pi + data.Pd);
    data.priors = d_read_from_file(argv[1], data.Nd);
    data.r = int_read_from_file(argv[3], data.Nr);
    // fprintf(stderr, "Got size of r %d \n", sizeof(&data.Nr) / sizeof(&data.Nr[0]));
    data.w = int_read_from_file(argv[5], data.N);
    data.d = int_read_from_file(argv[10], data.Nd);
    data.info = int_read_from_file(argv[4], data.Nd);

    // for (int i = 0; i < data.Nd; i++){
    //     fprintf(stderr, "Got value %d", data.info[i]);
    // }
    
    // long double h = backward(&data, 5, 5);
    //fprintf(stderr, "Got value %12.5e \n", h);
    
    decode(&data, llrf);
    // double t = logaddexp(-2, -INFINITY);
    // fprintf(stderr, "Test, got %12.5e \n", t);

    return 0;

}


