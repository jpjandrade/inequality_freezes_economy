#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>

double price(int c, int * goods_type);
void run_trading_instance(int N, double * capital, double * cash, double * pct_success, double ** ownership_average, double beta, double capital_goods_ratio, int dg_clever, gsl_rng * r, int M_user, int hist_flag, int thermalized);
double Mk_relative(int k);
void distribute_goods_random(int N, int M, int * goods, int * goods_type, double * cash);
void distribute_goods_clever(int N, int M, int M_K, int * goods, int * goods_type, double * cash);
void calc_goods_type(int * goods_type, int M, int M_K);
void generate_capital(int N, double * capital, double beta, gsl_rng * r);
void generate_adjusted_capital(int N, double * capital, double beta, gsl_rng * r);
void generate_bimodal_capital(int N, double * capital);
void generate_staircase_capital(int N, double * capital, double beta);
void read_from_state_files( double * capital, double * cash, int * goods, int * good_type, const int M, const int M_user, const int N);
int attempt_trade(int N, int M, int * goods, int * goods_type, double * cash, int ** Z, gsl_rng * r);
void calc_ownership(int N, int M, int * goods, int * goods_type, double * capital, int verbose, int verbose_file, float beta, float capital_goods_ratio);
void build_Z(int M, int N, int * goods, int * goods_type, int ** Z);
void update_ownership(int N, int ** Z, double ** ownership_average);
void update_histogram(int N, int M, int * goods, int * goods_type, int ** own_hist);
int compare (const void * a, const void * b);

int compare (const void * a, const void * b){
  if(*(const double*)a < *(const double*)b)
        return 1;
  return -1*(*(const double*)a > *(const double*)b);
}

