/**
 * simulation for Yegor Bazykin 2015 article
 * 
 * 2015-07-17
 * static landscape, tree structure events 000 ... 111
 * 
 * Fitness landscape is shaked at Poisson time
 * Substitutions occur at Poisson time
 * We detect events A->B->A, A->B>C, A->B->B at 0,1mln, 2mln,3mln,4mln,5mln,6mln.
 * on the static landscape the R/A statistic grows linearly with millions
 * from 1.60 - 1.70 no clear picture
 */

#include <iostream>
#include <random>
#include <vector>
#include <string>
#include <stack>
#include <sstream>
#include <chrono>
#include <set>
#include <algorithm>
#include <cmath>
#include <climits>

using namespace std;

bool logger=false;

string sAminoacids("GPAVLIMCFYWHKRQNEDST");
// it is a parameter we can change
// reversals shouls be more frequent with lesser states
int amino_acids_number=20;
//segment length is set automatically to allow 1 substitution per segment on average
int segment_length= 10000;
//there is no need for mu, because mu is in mutation time distribution
//double mu=0.000001;
//lesser population size leads to fixation probability increase
int N = 1000;
/**
 *     command line parameterers
 */
double substitution_delay_mean; //1000

//generator for normal distribution
unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine fitness_generator(seed1);

//main statistics
//stat for branches
int stat_step[6] = {0,0,0,0,0,0};
int stat_global_substs[3] = {0,0,0};
int stat_global_all[3] = {0,0,0};
int stat_segments[8]={0,0,0,0,0,0,0,0};

// a landscape for signle amino acid
class AAcid
{   
    public:
	AAcid(char c,double f);
	double fitness;	
	char A;
};

/**
 *  amino acid with a static landscape
 */
AAcid::AAcid(char c,double f)
{
     A=c;     
     fitness=f;
}

/**
*     Amino acids landscape in time and substitutions
*/
class Evolution
{
    vector<AAcid> landscape;    
    vector<int> times;
    
    vector<char>states;
    vector<int>state_numbers;
    
    //flat distribution to get random amino acid for substitution
    uniform_int_distribution <int> * i_dist;
    //flat real distribution for fixation
    uniform_real_distribution<double> * r_dist;	      
    
    public:
        vector<double>fitnesses;
	Evolution(double dfe_var);
	~Evolution();
	void evolve_tree();
	string sHistory;
	char get_state(int time);
	double  get_fitness(int time);
	int substitution();
	
	//distribution of fitness effects variation
	//when variance is small = 0.001 - nearly neutral reversals are frequent
	double dfe_var;
	//one distribution for all amino acids
	normal_distribution<double> * dfe;
	char current_state;
	double current_fitness;
};

Evolution::Evolution(double dfe_var)
{
 
    dfe = new normal_distribution<double>(1.0,dfe_var);  
    
    i_dist = new uniform_int_distribution <int> (0,amino_acids_number-1);
    r_dist = new uniform_real_distribution<double> (0.0,1.0);	    
  
    int i;
    //generate landscape
    if(logger)
	cout<<"========== Landscape =========="<<endl;
    for (i=0;i<amino_acids_number;i++)
    {  
	AAcid aa(sAminoacids[i],(*dfe)(fitness_generator));
	landscape.push_back(aa);
	if(logger)
	  cout<<aa.fitness<<endl;	
    }   
    
    //random amino acid for the first state
    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator1(seed1);      
    int iFirstState = (*i_dist)(generator1)%amino_acids_number;
    if(logger)
    {
      cout<<"========== Evolution =========="<<endl;
      cout<<"First state "<<iFirstState<<"\t"<<sAminoacids[iFirstState]<<endl;
    }
    current_state=landscape[iFirstState].A;    
    current_fitness=landscape[iFirstState].fitness;                
}

Evolution::~Evolution()
{
    delete dfe;    
}

int Evolution::substitution()
{
    int result = 0;
    unsigned seed2 = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator2(seed2);      
    
    int iNewState = (*i_dist)(generator2);
    double new_fitness=landscape[iNewState].fitness;
    
    //selection coefficient
    double s=new_fitness-current_fitness;
    //fixation classic formula 2.19 from Li. Molecular evolution
    double fixation_threshold = (1-exp(-2.0*s))/(1-exp(-4.0*N*s));
    //fixation probability ~ dfe_var;
    //segment length for 1 substitution / substitution_delay_mean * dfe = 1
    //segment length = substitution_delay_mean /2dfe ~ 1000 / *0.001 = 1 000 000
    double p = (*r_dist)(generator2);
    if (p<=fixation_threshold)//fix
    {
	if(logger)
	  cout<<"Fixed: "<<current_state<<"=>"<<landscape[iNewState].A<<"\t"<<s<<endl;
	current_state =  landscape[iNewState].A;
	current_fitness = landscape[iNewState].fitness;	  
	result=1;
    }	  
    return result;
}

void Evolution::evolve_tree()
{
    int iTimeFlow=0;
    int iTimeFlow1=0;
    int iTimeFlow2=0;
    int iTimeFlow3=0;
    int iTimeFlow4=0;
    int iTimeFlow5=0;
    int iTimeFlow6=0;
    
    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine substitution_generator(seed1);
    poisson_distribution<int> d_times(substitution_delay_mean);    
    
    int nsub=0;
    
    // segments 2,4,6
    //if(logger)
    //cout<<"========== Segment 2,4,6 =========="<<endl;
    /*
    while(iTimeFlow <= 3*segment_length)
    {       
	iTimeFlow += d_times(substitution_generator);	
	//cout<<iTimeFlow<<endl;
	nsub+=substitution();	
    }
    */
    nsub=0;
    
    if (!nsub) //no substitutions on segments 2-4-6
    {
	if(logger)
	    cout<<"========== Fake Segment =========="<<endl;
	char prev_state = current_state;
	char prev_fitness = current_fitness;
	iTimeFlow=0;
	//burning
	while(iTimeFlow <= 10*segment_length)
	{       
	    iTimeFlow += d_times(substitution_generator);	
	    nsub = substitution();	    
	}
      
	nsub=0;
	if(logger)
	    cout<<"========== Segment 1 =========="<<endl;
	current_state=prev_state;
	current_fitness=prev_fitness;    	
	iTimeFlow=0;
	while(iTimeFlow <= segment_length)
	{       
	    iTimeFlow += d_times(substitution_generator);	
	    nsub = substitution();
	    if (nsub)
		stat_step[0]=1;	
	}	
	
	current_state=prev_state;
	current_fitness=prev_fitness;    
	iTimeFlow=0;
	//segment 3
	if(logger)
	    cout<<"========== Segment 3 =========="<<endl;
	while(iTimeFlow <= segment_length)
	{       
	    iTimeFlow += d_times(substitution_generator);	
	    //cout<<iTimeFlow<<endl;
	    nsub=substitution();
	    if(nsub)
	      stat_step[2]=1;	
	}
	current_state=prev_state;
	current_fitness=prev_fitness;
	//segment 5
	if(logger)
	    cout<<"========== Segment 5 =========="<<endl;
	iTimeFlow=0;
	while(iTimeFlow <= segment_length)
	{       
	    iTimeFlow += d_times(substitution_generator);	
	    //cout<<iTimeFlow<<endl;
	    nsub=substitution();
	    if(nsub)
		stat_step[4]=1;	
	}
    
        int index = stat_step[4]*4+stat_step[2]*2+stat_step[0];
	//cout<<stat_step[4]<<"\t"<<stat_step[2]<<"\t"<<stat_step[0]<<endl;
	//cout<<index<<endl;
	stat_segments[index]++;	
    }
    
    for (int i=0;i<6;i++)
	stat_step[i]=0;
}

char Evolution::get_state(int time)
{
    char result;
    for (int i=0;i<times.size();i++)
    {
        if (times[i]>time)
	    break;
	else
	    result=states[i];
    }    
    return result;
}

double Evolution::get_fitness(int time)
{
    int result;
    for (int i=0;i<times.size();i++)
    {
        if (times[i]>time)
	  break;
	else
	  result=fitnesses[i];
    }    
    return result;
}

void print_results()
{
    for (int i=0;i<8;i++)
      cout<<stat_segments[i]<<endl;
}

int main(int argc, char *argv[])
{
    int replicates;
    double dDfeVar;
    if (argc==4)
    {
	basic_istringstream <char> ss1(argv[1]),ss2(argv[2]),ss3(argv[3]);
	ss1>>substitution_delay_mean;
	ss2>>dDfeVar;	
	ss3>>replicates;		
	segment_length = (int) (4.0*substitution_delay_mean/dDfeVar);
	//cout<<"Substitution_delay_mean = "<<substitution_delay_mean<<endl;
	//cout<<"Segment length = "<<segment_length<<endl;
    }
    else
    {
	cout<<"Usage: epistasis [substitution_delay_mean] [fitness_effects_var] [replicates]"<<endl;
	return 0;
    }  
    
    for (int i=0;i<replicates;i++)
    {
	Evolution ev(dDfeVar);    
	ev.evolve_tree();	
    }	
    print_results();	
    
    return 0;
}
