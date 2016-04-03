#include "gran_econ.h"

#define EPS 1e-7

//parameters
#define PI_0 0.01 //price of cheapest good (pi_k = PI_0 * G^k)
#define K 2 // number of different types of goods
#define G 2 //granularity of goods
#define C_MIN 1.    // minimal capital (not to bea changed)
#define alpha 1.    // the exponent of the price distributioin is alwyas 1

// variables used to produce staircase-like distributed Pareto distributions:
#define BASE 2
#define MathcalN 2
// alternatice C_MIN: #define C_MIN (PI_0 * pow(G,K-1)) //C_MIN is price of most expensive good

#define I_FOLLOW 0

//Dynamic related
#define T_TRIAL_STEPS 25000          // actual number of attempts at buying/selling that will be performed. (in millions). This is quite precisely proportional to the cpu time needed.

int main(int argc, char *argv[]){

  if (argc <6){
    printf("ERROR: need 5 arguments: N, beta, e, flag_thermalized, flag_hist  [M value of previous simulation]\n");
    return -1;
  }



  ///COMMAND LINE ARGS
  int N = atoi(argv[1]);
  double beta = atof(argv[2]);
  double capital_goods_ratio = atof(argv[3]); //capital_goods ratio ("e")

  int thermalized = atoi(argv[4]);
  int hist_flag = atoi(argv[5]);
  

  int M_user=-1; //sets the code to not load states from files...
  if(argc >6) M_user = atoi(argv[6]); //... unless we override it

  int dg_clever=1; //use block assignment of goods by default, 0 is random (slow).

  int i, k, c;
  double *cash, *capital, *pct_success;
  double **ownership_average;
  double ps_avg[K], ps_sq_avg[K];
  char out_name[200];
  FILE * out;
  clock_t start, end;
  start = clock();


  //initializes gsl_rng seed for the random functions (this is mostly irrelevant, one can ignore all this)
  const gsl_rng_type * T;  
  gsl_rng * r;
  double seed;
  gsl_rng_env_setup();  
  T = gsl_rng_default;
  r=gsl_rng_alloc (T);
  seed = (int) time (NULL) * getpid();
  //  seed = 42 ;   //if you want non random results uncomment this
  gsl_rng_set(r,seed);


  //memory alloc for stuff that depend on N and will be used inside run_trading_session
  cash = malloc(N*sizeof(double));
  capital = malloc(N*sizeof(double));

  pct_success = malloc(K*sizeof(double));
  ownership_average = malloc(N*sizeof(double*));
  for(i=0; i<N; i++) ownership_average[i] = malloc(K*sizeof(double));

  ///these are the observables that will be averaged over different capital distributions
  for(k=0; k<K; k++){ 
    ps_avg[k] = 0.;
    ps_sq_avg[k] = 0.;
  }


  printf("seed = %lf\n", seed);
  // generates initial capital, substitute the function as desired
  generate_adjusted_capital(N, capital, beta, r);    
  // generate_capital(N, capital, beta, r);  
  // generate_bimodal_capital(N, capital)
  // generate_staircase_capital(N, capital, beta);

  //sorts the capital, guy 0 is richest
  qsort(capital, N, sizeof(double), compare); 
  printf("Capital generated and sorted!\n");

  //MAIN FUNCTION CALL (the rest of the code is inside run_trading_instance)
  run_trading_instance(N, capital, cash, pct_success, ownership_average, beta, capital_goods_ratio, dg_clever, r, M_user, hist_flag, thermalized);

  //agregates observables
  for(k=0; k<K; k++){
    ps_avg[k]+=pct_success[k];
    ps_sq_avg[k] += pow(pct_success[k],2);
  }

  //output own_avg to file
  float GG = G;
  sprintf(out_name, "ownership_averages_K=%d_pi0=%.4f_g=%.1f_beta=%.2f_e=%.2f_N=%d_seed=%d.dat", K, PI_0, GG, beta, capital_goods_ratio, N, (int) seed);
  out = fopen(out_name, "w");
  for(i=0; i<N; i++){
    fprintf(out, "%f ", capital[i]);
    for(k=0; k<K; k++){
fprintf(out, "%f ", ownership_average[i][k]);
    }
    fprintf(out, "\n");
  }
  fclose(out);

  
  double total_capital = 0.0;
  for(i=0; i<N; i++){  
    total_capital += capital[i];
  }
    
  sprintf(out_name, "ps_K=%i_N=%i_beta=%.2f_e=%.2f.dat", K, N, beta, capital_goods_ratio);
  out = fopen(out_name, "a");
  for(k=0; k<K; k++) fprintf(out, "%.10f ", ps_avg[k]);
   fprintf(out, " %f %f %i %.5f %.1f %f \n", beta, capital_goods_ratio, N, PI_0, GG, (float) total_capital);
/*  beta=%.2f_e=%.2f.*/
  fclose(out);
  

  printf("--> Fraction of trades accepted (p_s): ");
  for(k=0; k<K; k++) printf("%f +- %f, ", ps_avg[k], sqrt(ps_sq_avg[k] - pow(ps_avg[k],2)));
  printf("\n");

  //maintenance
  free(capital);
  free(cash);
  free(pct_success);
  for(i=0; i<N; i++) free(ownership_average[i]);
  free(ownership_average);

  end = clock();
  printf("********Time elapsed: %.2fs\n", ((double)(end - start)) / CLOCKS_PER_SEC);
  return 0;
}

void run_trading_instance(int N, double * capital, double * cash, double * pct_success, double ** ownership_average, double beta, double capital_goods_ratio, int dg_clever, gsl_rng * r, int M_user, int hist_flag, int thermalized){
  

  int verbose_screen=0;
  int verbose_file=0;
  int T_MAX;
  int ownership_average_count=0;
  int i, j, l, k, c,t_M, T_Term;
  int z1, z2;
  int t;
  int good_traded;
  int total_trades, sum_trades, successful_trades[K], attempted_trades[K];
  double u, pi, total_capital, ps_inst;
  double * cum_capital;
  double ** followed_hist;
  int *goods, *goods_type;
  int ** Z;
  int hist_size[2]; 
  char out_name[400], ps_ts_name[400], in_name[200];
  FILE * out, *ps_time_series;
  int M, Mcontrol;
  //builds the ownership matrix
  Z = malloc(sizeof(int*) * N);
  for(i=0; i<N; i++) Z[i] = malloc(sizeof(int) * K);
  cum_capital = malloc(sizeof(int) * (N+1));

  double totalGoodsValue_perM;
  float GG = G;

  //we first assign cash = capital and calculate total capital
    if (M_user > -1){
        M=M_user;
    }else{
      total_capital = 0;
      cum_capital[0] = 0;
      for(i=0; i<N; i++){  
        cash[i] = capital[i];
        total_capital += capital[i];
	cum_capital[i+1] = total_capital;
      }
      assert(cum_capital[1] == capital[0]);
      assert(cum_capital[N] == total_capital);
      totalGoodsValue_perM = 0.;
      for(k=0; k<K; k++){
        totalGoodsValue_perM += Mk_relative(k) * PI_0 * pow(G, k);
      }
      Mcontrol = (int) ( total_capital / (capital_goods_ratio * totalGoodsValue_perM) );  
      //with Mcontrol and G we know all M_k and so we know M
      M = 0;
      for(k=0; k<K; k++){
        M+=(int) (Mcontrol*Mk_relative(k));
      }
      
    }
    printf("M defined = %d\n", M);
      /*goods is a vector of size M that, for each element c 
        will give you the owner goods[c] = 0, ..., N-1

    goods_type is a vector of size M that for each element c
    will give you the good type goods_type[c] = 0, ..., K-1 */
    goods = malloc(M*sizeof(int));
    goods_type = malloc(M*sizeof(int));


    if (M_user > -1){
    read_from_state_files( capital, cash, goods, goods_type, M, M_user, N);

      total_capital = 0;
      cum_capital[0] = 0;
      for(i=0; i<N; i++){  
        total_capital += capital[i];
	cum_capital[i+1] = total_capital;
      }
      	
      totalGoodsValue_perM=0.;
      for(k=0; k<K; k++){
        totalGoodsValue_perM += Mk_relative(k) * PI_0 * pow(G, k);
      }
      Mcontrol = (int) ( total_capital / (capital_goods_ratio * totalGoodsValue_perM) );  
      //with Mcontrol and G we know all M_k and so we know M
      M = 0;
      for(k=0; k<K; k++){
        M+=(int) (Mcontrol*Mk_relative(k));
      }

      printf("M=%d, M_user=%d \n ", M, M_user);
      assert(M==M_user);
    }
    printf("Capital / goods state (maybe) restored, M_user = %d\n", M_user);

  /*the code below will prepare the histograms for each consumer,
    each histogram goes from 0 to m_i and is calculated only for the *cheapest* good
  */
  if(hist_flag == 1){
    hist_size[0] = (int)(capital[I_FOLLOW]/(PI_0 * pow(G, 0))) + 1;
    if(K>1) hist_size[1] = (int)(capital[I_FOLLOW]/(PI_0 * pow(G, 1))) + 1;
    else hist_size[1] = 1;

    followed_hist = malloc(sizeof(double*) * hist_size[0]);
    for(z1=0; z1<hist_size[0]; z1++){
      followed_hist[z1] = malloc(sizeof(double) * hist_size[1]);
    }
  
    printf("Histogram allocated! Now we initialize it! Size: %d by %d\n", hist_size[0], hist_size[1]);
    for(z1=0; z1<hist_size[0]; z1++)
      for(z2=0; z2<hist_size[1]; z2++)
	followed_hist[z1][z2] = 0.;
    
    printf("Histogram of I_FOLLOW created\n");
  }
  
  if (M_user==-1){
  //create goods_type vector, goods_type[c] is the type of good c (0 for cheapest, K-1 for exp)
  calc_goods_type(goods_type, M, Mcontrol);
  printf("Now we will start distributing goods:\n");
  if(dg_clever==0) distribute_goods_random(N, M, goods, goods_type, cash);
  if(dg_clever==1) distribute_goods_clever(N, M, Mcontrol, goods, goods_type, cash);
  }

  printf("Building Z\n");
  build_Z(M, N, goods, goods_type, Z);
  //calculates and prints initial ownership on screen (1st binary argument) or on file (2nd)
  //calc_ownership(N, M, goods, goods_type, capital, verbose_screen, verbose_file, beta, capital_goods_ratio);

  sprintf(ps_ts_name, "ps_time_series_K=%d_pi0=%.4f_g=%.1f_beta=%.2f_e=%.2f_N=%d_M=%d.dat", K, PI_0, GG, beta, capital_goods_ratio, N, M);
  ps_time_series = fopen(ps_ts_name, "w");
  //debug statement
  printf("---INITIAL CONFIGURATION---\n");
  printf("M = %d,", M);
  for(k=0; k<K; k++) printf(" M_%d = %d", k, (int) (Mcontrol*Mk_relative(k)) );
  printf(", N = %d\n", N);
  printf("Total capital: %.2f\n", total_capital);
  printf("Total good value: %.2f\n", totalGoodsValue_perM * Mcontrol );


  //initialize variables
  for(i=0; i<N; i++){
    for(k=0; k<K; k++){
      ownership_average[i][k]=0.;
    }
  }

  T_MAX = (int) ((long) T_TRIAL_STEPS * (long) 1000000 / (long) M) ;

  for(k=0; k<K; k++) pct_success[k] = 0;
  if(thermalized) T_Term = (int)(0.2*T_MAX);
  else T_Term = (int)(0.8 * T_MAX); //just get the last 20% points
  
  /////MAIN LOOP OF THIS FUNCTION
  //each t is a Monte Carlo step, runs over M*N random goods
  for(t=0; t<T_MAX; t++){
    if((10*t)%T_MAX==0) printf("%.0f%% done\n", (double)(100*t)/(T_MAX));
    
    //TRADING SESSIONS
      //we do observable updating inside this loop to avoid overflow
      total_trades = 0; 
      for(k=0; k<K; k++) successful_trades[k] = 0; 
      for(k=0; k<K; k++) attempted_trades[k] = 0;

      for(t_M = 0; t_M < M; t_M++){ //separate loop to avoid overflow and (mostly) to do some averages on the run.
	
	    /*here is where the actual trades happen: 
	      good_traded will return c from 0 to M-1 in case good c is chosen and traded
	      it will return c from -1 to -M in case good -c-1 (sorry) is chosen 
	      and fails to trade */
	good_traded = attempt_trade(N, M, goods, goods_type, cum_capital, cash, Z, r);
	    
	
	    if(t == T_Term && t_M ==0){printf("starting the recording - as if we were thermalized: now.\n");}
	    if(t >= T_Term){ //data collection, after termalization

	      /*this updates the ownership_average statistic and histograms
	        at each M time step to avoid autocorrelation */
	      
	      //if( (t_M % SAMPLE_STEPS ==SAMPLE_STEPS-1) || (t==T_MAX-2 && t_M == 0 ) ){  // WARNING
	      if(t_M == 0){
	        //upgrades avg ownership
	        update_ownership(N, Z, ownership_average);
	        ownership_average_count +=1;
	        //printf("%d %d\n", z_follow[0], z_follow[1]);
	        //upgrades instantaneous histogram of poorest type of goods
	        if(hist_flag == 1){
		        if(K>1) followed_hist[Z[I_FOLLOW][0]][Z[I_FOLLOW][1]]+=1e-6;
		        else followed_hist[Z[I_FOLLOW][0]][0] += 1e-6;
		        }
	      }
       
	      //here we update for every trade to calculate p_s and other obs
	      if(good_traded >=0){
	        successful_trades[goods_type[good_traded]]++;
	        attempted_trades[goods_type[good_traded]]++;

	      }
	      else{ //trade failed, we have to do c = -(output of function + 1) to find good type
  
	        attempted_trades[goods_type[-good_traded -1]]++;
    
	      }
        
	    } //close if (t>t_term)
      } //closes t_M loop

      //here we collect the trade totals to fractions so we avoid overflows
      sum_trades = 0;
      for(k=0; k<K; k++) sum_trades += attempted_trades[k];
      assert(sum_trades == 0 || sum_trades == M); //make sure we didn't forget to count anything


      if(sum_trades > 0){
	    for(k=0; k<K; k++){
	      if(attempted_trades[k]>0){
	        ps_inst = successful_trades[k]/((double)attempted_trades[k]);
	        fprintf(ps_time_series, "%f ", ps_inst);
	        pct_success[k] += ps_inst;
	      }
	      else{
	        pct_success[k] += 1.;
	        printf("Trade skipped, be careful, numbers may be skewed if happens often!\n");
	      }
	    }
      fprintf(ps_time_series, "\n ");
      }

  }

  for(k=0; k<K; k++) pct_success[k] /= (T_MAX - T_Term) ; 
  for(i=0; i<N; i++){
    for(k=0; k<K; k++){
      ownership_average[i][k]/= (ownership_average_count);
    }
  }

  printf("recording the system parameters and state \n");
  sprintf(out_name, "capital_M=%d", M);
  out = fopen(out_name, "w");
  for(i=0; i<N; i++){
    fprintf(out, "%.10f %.10f\n", capital[i], cash[i] );
  }
  fclose(out);
  sprintf(out_name, "goods_M=%d", M);
  out = fopen(out_name, "w");
  for(i=0; i<M; i++){
    fprintf(out, "%d %d\n", goods[i], goods_type[i] );
  }
  fclose(out);
  printf("recording completed. \n");
  
  //outputs histograms to file, the rest will be collected in main function
  if(hist_flag == 1){
    sprintf(out_name, "histogram_K=%d_pi0=%.4f_g=%.1f_beta=%.2f_e=%.2f_N=%d_M=%d_agent=%d.dat", K, PI_0, GG, beta, capital_goods_ratio, N, M, I_FOLLOW);
    out = fopen(out_name, "w");
    for(z1=0; z1< hist_size[0]; z1++){
      if(K>1){
	for(z2=0; z2 < hist_size[1]; z2++){	  
	  fprintf(out, "%d %d %e \n", z1, z2, followed_hist[z1][z2]);
	}
      }
      else
	fprintf(out, "%d %e \n", z1, followed_hist[z1][0]);
    }
    fclose(out);
  }

  
  fclose(ps_time_series);
  

  
  
  printf("going to clear some stuff here\n");
  //maintenance
  free(goods);
  free(goods_type);
  if(hist_flag == 1){
    for(z1=0; z1<hist_size[0]; z1++){
      free(followed_hist[z1]);
    }
    printf("freed followed hist\n");
    free(followed_hist);
    
    printf("freed followed hist2\n");
  }

  for(i=0; i<N; i++) free(Z[i]);
  free(Z);
}

void generate_capital(int N, double * capital, double beta, gsl_rng * r){
  int i,k;
  for(i=0; i<N; i++){
    capital[i] = gsl_ran_pareto(r, beta, C_MIN);
  }
}

void generate_adjusted_capital(int N, double * capital, double beta, gsl_rng * r){
  int i,k;
  for(i=0; i<N; i++){
    capital[i] = gsl_ran_pareto(r, beta, C_MIN);
  }
  int agent;
  double totCap;
  printf("capital has to be adjusted\n");
  totCap=0.; for(i=0; i<N; i++){ totCap += capital[i]; }
  if( totCap/N < (beta/(beta-1.)) ) {
        capital[0] += ((beta/(beta-1.)) - totCap/N)*N;
        printf("capital has been increased\n");
  }else{
      agent = 0 ;
      qsort(capital, N, sizeof(double), compare); 
      while ( totCap/N > (beta/(beta-1.)) ) {
         capital[agent] -= ( totCap/N - (beta/(beta-1.)) )*N;
         if ( capital[agent] < C_MIN ){
            capital[agent]=1.;
            agent+=1;
         }
         totCap=0.; for(i=0; i<N; i++){ totCap += capital[i]; }
     }
     printf("capital has been decreased (rare event), we had to go up to agent #%d\n", agent);
  }
}


void generate_bimodal_capital(int N, double * capital){
  int i,k;
  for(i=0; i<N/2; i++){
    capital[i] = 1.00001;
  }  
  for(i=N/2; i<N; i++){
    capital[i] = 3.00001;
  }  
}


void generate_staircase_capital(int N, double * capital, double beta){
    int I,i, NNi,j, Ninit,flagExit=0;
    Ninit=N;
    float cci;
    I=0;
    while (I<N-1) {
        for(i=0; i<MathcalN; i++){
            cci = pow(BASE, i+0.5);
            NNi = (int)( Ninit * pow(BASE, -i*beta) * (1.0 - pow(BASE, -beta) ));
/*            printf("I=%d   i=%d   NNi=%d    cci=%f \n", I, i, NNi,  cci);*/
            for(j=I; j<I+NNi; j++){
                capital[j] = cci;
            }
            I+=NNi;
            if(NNi<10){flagExit=1;}
            if(I>N-10){I=N-1; break;}
        }
        Ninit = N-I;
        if (flagExit==1){break;}
    }
    printf("I=%d   i=%d   - there are nn=N-I=%d  sites that are excluded \n", I,i, N-I);
    for(j=I; j<N; j++){
        capital[j] = 0.0; // we eliminate these players, to keep the ratios correct.
    }
}


void read_from_state_files( double * capital, double * cash, int * goods, int * goods_type, const int M, const int M_user, const int N){

    int i;
    char in_name[100];

    double a;
    float * tempooo;
    float * trash;
    trash = malloc(2*sizeof(float));
    tempooo = malloc(2*sizeof(float));
    FILE *fp;
    if (M_user > -1){
        sprintf(in_name, "capital_M=%d", M_user);
        fp = fopen(in_name,"r"); // read mode
        if( fp == NULL )  { perror("Error while opening the file 'capital'.\n"); exit(EXIT_FAILURE); }
        for (i = 0; i < N; i++)
        {
            a=fscanf(fp, "%10f ", &tempooo[0]);
            a=fscanf(fp, "%10f ", &trash[0]);
            capital[i] = (double) tempooo[0];
            a=fscanf(fp, "%10f ", &tempooo[1]);
            cash[i] = (double) tempooo[1];
            a=fscanf(fp, "%10f ", &trash[1]);
        }
        fclose(fp);
	sprintf(in_name, "goods_M=%d", M_user);
        fp = fopen(in_name,"r"); // read mode
        if( fp == NULL )  { perror("Error while opening the file 'goods'.\n"); exit(EXIT_FAILURE); }
        for (i = 0; i < M; i++)
        {
            a=fscanf(fp, "%10f ", &tempooo[0]);
          //  a=fscanf(fp, "%10f ", &trash[0]);
            goods[i] = (int) tempooo[0];
            a=fscanf(fp, "%10f ", &tempooo[1]);
            goods_type[i] = (int) tempooo[1];
          //  a=fscanf(fp, "%10f ", &trash[1]);
        }
        fclose(fp);
    }
    free(tempooo);
    free(trash); 
}


void distribute_goods_random(int N, int M, int * goods, int * goods_type, double * cash){
  int i,c;
  double pi;
  int allowed_buyer;
  for(c=M-1; c>=0; c--){ //we start from most expensive goods
    goods[c] = -1;
    pi = price(c, goods_type);
    //printf("giving out good %d, price: %f\n", c, pi);
    allowed_buyer =N-1;
    for(i=allowed_buyer; i>=0; i--){

      //printf("dude %d, cash = %f: ", i, cash[i]);
      if(cash[i] >= pi){
	cash[i] = cash[i] - pi;
	goods[c] = i;
	//printf("SOLD! New cash: %f\n", cash[i]);
	/* if(i == i_follow){ */
	/*   kind = goods_type[c]; */
	/*   z_follow[kind] ++; */
	/* } */
	i=-1;
      }else{
	allowed_buyer --;
      }
    }
  }
}

void distribute_goods_clever(int N, int M, int Mcontrol, int * goods, int * goods_type, double * cash){
  int i,c,d, k;
  int next_flag;
  int verbose = 0;
  int M_k, m_i, remaining_goods;
  double pi;
  //We start with the poorest guy:
  i = N-1;

  //We start with the cheapest good
  c = 0;
  for(k=0; k<K; k++){
    if(k>0) i=0;
    pi = PI_0 * pow(G,k);
    M_k = (int) (Mcontrol*Mk_relative(k)); 
    remaining_goods = M_k;

    while(remaining_goods>0){
      next_flag = 0;

      m_i = floor(cash[i]/pi);
      //printf("We have %d goods of type %d to still give. Guy %d is on, with cash %f and m %d\n", remaining_goods, k, i, cash[i], m_i);
      if(m_i < remaining_goods){ //buys m_i goods
	    next_flag = 1;	
      }
      else{
    	m_i = remaining_goods;
      }
      assert(m_i>=0);
      cash[i]-= m_i * pi;
      for(d=c; d<c+m_i; d++) 
      {
        goods[d] = i;
      }
      c = c+m_i;
      remaining_goods -= m_i;
      if(verbose) printf("Guy %d just bought %d goods of type %d! He has cash %f and we still have %d goods to go!\n", i, m_i, k, cash[i], remaining_goods);
      if(verbose) printf("c=%d, m_i = %d\n", c, m_i);
      assert(cash[i]>=0);
      assert(remaining_goods >= 0);
      if(m_i>0) assert(goods_type[c-1] == k);
      if(next_flag){
	if(k == 0) i--;
	else i++;
      }
      assert(i>=0);
    }
  }
}



int attempt_trade(int N, int M, int * goods, int * goods_type, double * cum_capital, double * cash, int ** Z, gsl_rng * r){
  
  int i, j, c, k, s, kind;
  int i_found, print_flag = 0;
  double pi;
  double u, cap;
  c = gsl_rng_uniform_int(r,M);
  j = goods[c]; //owner of good c
  do{
    // i = gsl_rng_uniform_int(r, N);  //old, random choice
    k = gsl_rng_uniform_int(r, 50) + 1;
    s = gsl_rng_uniform_int(r, 2);
    s = 2*s - 1;
    assert(s == 1 || s == -1);
    i = j + s * k;
    if(i < 0){
      if(j==0) i = 1;
      else i = 0;
    }
    if(i > N - 1){
      if(j == N - 1) i = N - 2;
      else i = N - 1;
    }
    /* u = gsl_rng_uniform(r); */
    /* cap = u * cum_capital[N-1]; */
    /* i_found = 0; */
    /* i = 0; */
    /* while(i_found == 0){ */
    /*   if(cum_capital[i+1] < cap) i++; */
    /*   else i_found = 1;       */
    /* } */
    //printf("Agent chosen for trade: %d, cum_capital: %.3f - %.3f, cap = %.3f\n", i, cum_capital[i], cum_capital[i+1], cap);
  }while(i==j);
  pi = price(c, goods_type);
  if(cash[i] >= pi){
    cash[j] = cash[j] + pi;
    cash[i] = cash[i] - pi;
    goods[c] = i;
    
    //if(measure == 1){
        kind = goods_type[c];
        Z[i][kind]++;
        Z[j][kind]--;
    //}
    
    return c;
  }
  else{
    return -c-1; // -1 if good 0 failed to trade, -2 if good 1 failed to trade, etc
  }
}

void calc_ownership(int N, int M, int * goods, int * goods_type, double * capital, int verbose, int verbose_file, float beta, float capital_goods_ratio){
  int i, c, k;
  int own[K];
  char out_name[200];
  FILE * out;    
  float GG = G;
  if(verbose_file){
    sprintf(out_name, "ownership_K=%d_pi0=%.4f_g=%.1f_beta=%.2f_e=%.2f_N=%d_M=%d.dat", K, PI_0, GG, beta, capital_goods_ratio, N, M);
    out = fopen(out_name, "w");
  }
  if(verbose) printf("ownership: \n");
  for(i=0; i<N; i++){
    for(k=0; k<K; k++) own[k] = 0;
    
    if(verbose) printf("Consumer %d, init. capital %.2f: ", i, capital[i]);
    for(c=0; c<M; c++){
      if(goods[c]==i){
	own[goods_type[c]]++;
      }
    }
    if(verbose){
      for(k=0; k<K; k++) printf("%d ", own[k]);
      printf("\n");
    }
    if(verbose_file){
      fprintf(out, "%f ", capital[i]);
      for(k=0; k<K; k++) fprintf(out, "%d ", own[k]);
      fprintf(out, "\n");
    }
    
  }  
  if(verbose_file) fclose(out);
}

void build_Z(int M, int N, int * goods, int * goods_type, int ** Z){
  int i, j,c, l, k;
  for(i=0; i<N; i++){
    for(k=0; k<K; k++){
      Z[i][k]=0;
    }
  }
  
  for(c=0; c<M; c++){
    j = goods[c];
    l = goods_type[c];
    Z[j][l]++;
  }
}

void update_ownership(int N, int ** Z, double ** ownership_average){
  int i,k;
  for(i=0; i<N; i++){
    for(k=0; k<K; k++){
      ownership_average[i][k]+=Z[i][k];
    }
  }
}


void calc_goods_type(int * goods_type, int M, int Mcontrol){
  int c, k, res_type;
  int type_vec[K];
  for(k=0; k<K; k++){
    type_vec[k] = (int) (Mcontrol*Mk_relative(k));
    if(k>0) type_vec[k] += type_vec[k-1];
  }
  assert(type_vec[K-1] == M);

  for(c=0; c<M; c++){
    res_type = -1;
    k=0;
    while(res_type == -1 && k<K){
      if(c < type_vec[k])
	res_type = k;
      else
	k++;
    }
    goods_type[c] = res_type;
  }      
}
 

double Mk_relative(int k){
    return (1.0-pow(1.0/G,alpha)) * pow(1.0/G, k*alpha) ; 
}


double price(int c, int * goods_type){
  int k_type, good_type;
  double pi;
  k_type = goods_type[c];
  pi = PI_0 * pow(G,k_type);
  return pi;
}


  /* sprintf(out_name, "expensive_traded_ec%d.dat", (int)(ec_mult*10)); */
  /* out = fopen(out_name, "w"); */
  /* fprintf(out, "%f %f\n", extra_cash, fr_expensive_traded/((double)M_2 / M)); */
  /* fclose(out); */
