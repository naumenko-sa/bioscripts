/* 10 March 2011 - 12 June 2011. Simulated evolution, for comparing MK and BK tests */

/*    This function simulates evolution of one diallelic (0 and 1) locus.

Using this function, we consider loci of 5 kinds: neutral (synonymous and nonsynonymous, but who cares)
under constant selection (say, Ns = 1, 10, and 100) and under changing selection (perhaps, Ns = 100).

The locus is characterized by coefficient of selection and waiting time before this coefficient switches sign.
Also, there is a mutation rate and actual (= effective, for the population is a WF one) population size.
Also, there is the phylogenetic tree of a prescribed topology.
Also, I pass a random seed into the function.
Output consists number of alleles 1 present currently in each species, and a randomly sampled allele from each species.

User: take care to run this function with different random seeds!

User: To determine initial probability of the allele, you may use equilibrium distribution.
Alternatively, start from the currently fittest one. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define MM 0.0000001    /* mutation rate - it more or less determines other parameters, as Nm must be 0.01, Ns must be reasonable, and the number of generations must produce the desired Ks */
#define NGEN 70000000L  /* number of generations */
#define NN 100000L      /* population size */
#define PMAX1 100L      /* maximal number of rare allele carriers, for which Poisson distribution is used */
#define PMAX2 300L      /* maximal number of rare allele carriers, produced by drift, modelled by Poisson distribution */
#define GMAX 1500L      /* must exceed 5(NN^0.5) - maximal change in the number of allele carriers, produced by drift, modelled by Gaussian distribution */

static void init_genrand(unsigned long);                 /* function 1 of random number generator */
static double genrand_real3(void);                       /* function 2 of random number generator */
static void poissoncr(FILE *, long, long double *);      /* making a Poisson distribution */
static void gausscr(FILE *, long double, long double *); /* making a half-Gaussian distribution */
static long LocusHistoryInit(FILE *, long double (*)[PMAX2+1], long double (*)[GMAX+1]);
static void LocusHistory(FILE *, FILE *, long double (*)[PMAX2+1], long double (*)[GMAX+1], long, long double, long, long *, long *, long *, long *); /* MY FUNCTION */

int main(void) /* This main program calls my function. It provides it with random seed. */
{
  long i, k, seed, switches, node[6], alleles[7], AllNum[7], a162[7]; /* There are 6 species in our phylogeny */
  long double s, R, (*poissonTrans)[PMAX2+1], (*gaussTrans)[GMAX+1];

  FILE *fout;  if ((fout  = fopen("/dev/null","w"))  == NULL) goto quit;
  FILE *fout1; if ((fout1 = fopen("PositiveSelectionSimulationF07.out","w")) == NULL) goto quit;
  if ((poissonTrans = calloc((PMAX1+1)*(PMAX2+1),   sizeof(long double))) == NULL) { fprintf(fout," No memory-1! "); goto quit; } /* Memory allocation */
  if ((gaussTrans   = calloc((NN-PMAX1+1)*(GMAX+1), sizeof(long double))) == NULL) { fprintf(fout," No memory-2! "); goto quit; } /* Memory allocation */

  printf (" zapolnyayu figushki!\n");
  i = LocusHistoryInit(fout,poissonTrans,gaussTrans); if (i == -1) goto quit;
  printf (" zapolnila figushki!\n");

/* times of internal nodes - reflect Drosophila tree - topology is defined inside my function. */
  node[0] = -53800000; /* anopheles */
  node[1] = -12100000; /* virilis */
  node[2] = -10700000; /* pseudoobscura */
  node[3] =  -7200000; /* ananassae */
  node[4] =  -1400000; /* erecta */
  node[5] =   -700000; /* simulans */


  time_t seconds;
  seconds = time (NULL); 

  seed = seconds;                    /* Initializing random number generator - each time, use a different seed */
  init_genrand(seed);



for (i = 0; i != 100; i ++)   /* calling my function many times */
{
  for (k = 0; k != 1; k ++)
  {
/* coefficient(s) of selection for the derived allele and waiting time for switches of the sign of coefficient of selection; negaive means no switches */
    if (k == 0) s = -0.01; else s = -0.1;      /* different coefficients of selection */
    if (k <= 3) switches = -1; else switches = 200000000;                                                   /* different switching of selection regimes */     

    printf  (       " %12ld %10ld %11.9f %5ld %10ld %10.8f %10ld\n", NGEN, NN, MM, i, seed, (float)s, switches);
    fprintf (fout,  " %12ld %10ld %11.9f %5ld %10ld %10.8f %10ld\n", NGEN, NN, MM, i, seed, (float)s, switches);
    fprintf (fout1, "%12ld\t%10ld\t%11.9f\t%5ld\t%10ld\t%10.8f\t%10ld\t", NGEN, NN, MM, i, seed, (float)s, switches);
 
/*  inputs:  two transition matrices, random generator seed, selection coefficient, waiting time for its switches, times of branchings.
    outputs: randomly chosen allele in each species, numbers of allele carriers in each species. */
    LocusHistory(fout, fout1, poissonTrans, gaussTrans, seed, s, switches, node, alleles, AllNum, a162);

    fprintf (fout1, "\t%7ld\t%7ld\t%7ld\t%7ld\t%7ld\t%7ld\t%7ld\t", AllNum[0], AllNum[1], AllNum[2], AllNum[3], AllNum[4], AllNum[5], AllNum[6]);
    fprintf (fout1, "%7ld\t%7ld\t%7ld\t%7ld\t%7ld\t%7ld\t%7ld\t", alleles[0], alleles[1], alleles[2], alleles[3], alleles[4], alleles[5], alleles[6]);
    fprintf (fout1, "%7ld\t%7ld\t%7ld\t%7ld\t%7ld\t%7ld\t%7ld\n", a162[0], a162[1], a162[2], a162[3], a162[4], a162[5], a162[6]);
  }
}

quit: fclose(fout); fclose(fout1); free(poissonTrans); free(gaussTrans);
  return 0;
}

/* Filling 2 transition matrices - this takes a lot of time and thus is done in a separate function */
long LocusHistoryInit(FILE *fout, long double (*poissonTrans)[PMAX2+1], long double (*gaussTrans)[GMAX+1])
{
  long LinCount, i, j;
  long double var, *poisson, *gauss;

  if ((poisson = calloc(PMAX2+1, sizeof(long double))) == NULL) { fprintf(fout," No memory-3! "); return -1; };
  if ((gauss   = calloc(GMAX+1,  sizeof(long double))) == NULL) { fprintf(fout," No memory-4! "); return -1; };

/*    Calculate transition matrices. These matrices describe only drift - selection, if any, comes later. */
/* Generating an array of Poisson distributions, for dynamics when an allele is very rare or very common */
  for (i = 1; i != PMAX1+1; i ++)
  {
    poissoncr(fout,i,poisson);                                        /* creating a Poisson distribution with parameter i */
    for (j = 0; j != PMAX2+1; j ++) poissonTrans[i][j] = poisson[j];  /* array of Poisson distributions with different parameters */
  }

/* Generating an array of half-Gaussian distributions, for dynamics when an allele has an intermediate frequency */
  for (i = PMAX1+1; i != NN-PMAX1; i++)
  {
    var = (float)i*(NN-i)/NN;                                     /* variance of a Gaussian distribution which corresponds to a particular i */
    gausscr(fout,var,gauss);                                      /* creating a half-Gaussian distribution with variance var */
    for (j = 0; j != GMAX+1; j ++) gaussTrans[i][j] = gauss[j];   /* array of half-Gaussian distributions with different variances */
  }
return 0;
}

/* My funcion works in the following way:
1) Loop by generations.
  a) Loop by lineages.
  b) Mutation - the number of derived alleles can increase and/or decrease by 1.
  c) If a lineage is monomorpic, fuck it.
  d) If a lineage is polymorphic, do drift and selection.
2) If this is the time of a cladogenesis, perform it, and increment the counter of lineages.
3) When present time is reached, report the results. 162? */
void LocusHistory(FILE *fout, FILE *fout1, long double (*poissonTrans)[PMAX2+1], long double (*gaussTrans)[GMAX+1], long seed, long double sinit, long switches, long *node, long *allele, long *AllNum, long *a162)
{
  long LinCount, i, j, t, itmp;
  long double R, var, s[7], tmp;

  AllNum[0] = 0; /* The root allele AllNum may be initiated stochastically, depending on s */
  s[0] = sinit;  /* Initial coefficient of selection */
  LinCount = 1;  /* Only one lineage, initially */

/*            REAL WORK          */
  for (t = -NGEN; t != 1; t ++) /* loop by generations */
  { 
    for (i = 0; i != LinCount; i ++) /* loop by lineages */
    {    
      R = genrand_real3();
      if (switches > 0 && R < 1.0/switches) /* switching the direction of selection */
      {
        s[i] = -s[i];
        fprintf (fout1, " (%5ld %10ld %10.8f) ", i, t, (float)s[i]);
      }

/* Changes of allele number due to mutation */
      R = genrand_real3(); if (R < MM*(NN-AllNum[i])) AllNum[i] ++; 
      R = genrand_real3(); if (R < MM*AllNum[i])      AllNum[i] --;

      if (AllNum[i] == 0 || AllNum[i] == NN) continue; /* if something is fixed, this is it for the lineage */

/* Changes of allele number due to drift */
      R = genrand_real3();
      if (AllNum[i] <= 100)            /* on lower edge do Poisson of rare allele */
      {
        for (j = 0; j != PMAX2+1; j ++) { if (R < poissonTrans[AllNum[i]][j]) break; }   /* find the Poisson number */
        AllNum[i] = j;      /* assign this number */
      }
      else if (AllNum[i] >= NN-100)    /* on higher edge do Poisson of rare (the other!) allele */
      {            
        for (j = 0; j != PMAX2+1; j ++) { if (R < poissonTrans[NN-AllNum[i]][j]) break; } /* find the Poisson number */
        AllNum[i] = NN-j;   /* assign this number */
      }
      else                             /* otherwise, do halfGaussian increment, with random choice of direction */
      {
        for (j = 0; j != GMAX+1; j ++) { if (R < gaussTrans[AllNum[i]][j]) break; }       /* find the Gaussian increment */                                           
        R = genrand_real3(); if (R > 0.5) AllNum[i] += j; else AllNum[i] -= j; /* use this increment */
      }

/* Changes of allele number due to selection - multiply by 1+s, do rounding using a random addition */
      tmp = (AllNum[i]*(1.0+s[i])/(AllNum[i]*(1.0+s[i])+NN-AllNum[i]))*NN;
      itmp = tmp; R = genrand_real3(); if (R < tmp-itmp) itmp ++;
      AllNum[i] = itmp;
    }

/* Internal nodes - in all cases, only one of the new lineages will bifurcate later */
    if      (t == node[0]) { AllNum[LinCount] = AllNum[LinCount-1]; s[LinCount] = s[LinCount-1]; LinCount ++; } /* anopheles branches off */
    else if (t == node[1]) { AllNum[LinCount] = AllNum[LinCount-1]; s[LinCount] = s[LinCount-1]; LinCount ++; } /* virilis branches off */
    else if (t == node[2]) { AllNum[LinCount] = AllNum[LinCount-1]; s[LinCount] = s[LinCount-1]; LinCount ++; } /* pseudoobscura branches off */
    else if (t == node[3]) { AllNum[LinCount] = AllNum[LinCount-1]; s[LinCount] = s[LinCount-1]; LinCount ++; } /* ananassae branches off */
    else if (t == node[4]) { AllNum[LinCount] = AllNum[LinCount-1]; s[LinCount] = s[LinCount-1]; LinCount ++; } /* erecta branches off */
    else if (t == node[5]) { AllNum[LinCount] = AllNum[LinCount-1]; s[LinCount] = s[LinCount-1]; LinCount ++; } /* simulans branches off */

    if (-t % 100000 == 0)   /* reporting progress - remove later */
    {
      printf (" %12ld ", t);
      for (i = 0; i != LinCount; i ++) printf ("%7ld ", AllNum[i]);
      printf ("\n");

      fprintf (fout, " %12ld ", t);
      for (i = 0; i != LinCount; i ++) fprintf (fout, "%7ld ", AllNum[i]);
      fprintf (fout, "\n");
    }
  }

/* After-run sampling of one allele from each species. What about 162? */
  for (i = 0; i != LinCount; i ++)
  {
    R = genrand_real3(); 
    if (R < (float)AllNum[i]/NN) allele[i] = 1; else allele[i] = 0;
    a162[i] = 0;
    for (j = 0; j != 162; j++) {R = genrand_real3(); if (R < (float)AllNum[i]/NN) a162[i] ++; }
  }
} /* END OF MY FUNCTION */

/* Generating Poisson distribution */
void poissoncr (FILE *fout, long Lambda, long double *poisson)
{
  long i;

  poisson[0] = 1.0;
  for (i = 1; i != PMAX2+1; i++) poisson[i] = poisson[i-1]*Lambda/i;
  for (i = 0; i != PMAX2+1; i++) poisson[i] *= exp(-Lambda);
  if (poisson[PMAX2] > 0.00000001) printf("poissoncr - Kozel1!");

  for (i = 1; i != PMAX2+1; i++) poisson[i] = poisson[i]+poisson[i-1];
  if (poisson[PMAX2] < 0.99999999) printf("poissoncr - Kozel2!");
  poisson[PMAX2] = 1.0;
}

/* Generating half-Gaussian distribution */
void gausscr (FILE *fout, long double var, long double *gauss)
{
  long i;
  long double Summ = 0.0;
  for (i = 0; i != GMAX+1; i++) { gauss[i] = exp(-(i*i)/(2*var)); Summ += gauss[i]; }
  for (i = 0; i != GMAX+1; i++) gauss[i] /= Summ;
  for (i = 1; i != GMAX+1; i++) gauss[i] = gauss[i]+gauss[i-1];
  if (gauss[GMAX] < 0.99999999) printf("gausscr - Kozel1!");
  gauss[GMAX] = 1.0;
}


/* BELOW, ARE TWO STANDARD FUNCTIONS */
/* Period parameters http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html - Mersenne twister */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}
