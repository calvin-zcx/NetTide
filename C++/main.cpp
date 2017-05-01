#include <stdio.h>
#include <string.h>
#include <map>
#include <hash_map>
#include <set>
#include <vector>
#include <time.h>
#include <algorithm>
#include <string>
#include <math.h>
#include <random>
#include <queue>
using namespace std;

/*
Author: Chengxi Zang
Date: 2017-05-01
Version: 2.0

Goal: Generator realistic growth dynamics, both for node and link, of social networks.

The growth dynamics of node and link are captured by NETTIDE equations:
	*Node adoption:
		n-dot = beta / t^theta * n * (N - n )
	*Link build between infected-infected and infected-newly infected:
		e-dot = beta_prime / t^theta * n * ( alpha * (n-1)^gamma - e/n ) + 2*n-dot
	*Assumption:
		beta and beta_prime : growth rate of nodes and links;
		1/t^theta : fizzling effect
		the alpha and gamma describe the average number of friends in the ego network.

Two generators:
	*NetTide-Process: captures the micro-level stochastic interactions within a network.
	*NetTide-Survival: captures the growth dynamics from hazard rates (hazard process).

Input 
	*Input network, reads the edge file 
	    Data format: blank-separated pairs , e.g 0 1 represents 0->1
		Network: e.g.,
			1.Complete network, random network
			2.Kronecker network
			3.Or any other network
        Nodes have been numbered 0-$n-1$
	*Parameters: beta, theta, beta_prime, alpha, gamma, N

Output:
    *n(t) over time 
	*e(t) over time. 
	*For NetTide-Process, detailed evovling process of nodes and links are recorded.
*/



// Global data structure
struct stStateOfInfection
{
	int userId;
	bool isInfected;
	int infectedByWho;
	float timestamp;
};

// Random number generator for decision of the infection.
default_random_engine generator(time(NULL));
uniform_real_distribution<double> distribution(0.0, 1.0);
// usage : double number; number = distribution(generator);


////PART1: Auxiliary functions 
void loadNetwork(int number_of_node, vector<vector<int>> &G, char *fileName)
{
/*
	G is either vUser_vFriend or vUser_vFriendForLink
	Let the node id from 0:number_of_node-1
	G[number_of_node][VARY]
*/
	FILE *fi; //fi = fopen("kronecker_3125.txt", "rb");
	fi = fopen(fileName, "rb");
	int p, c; 
	G.resize(number_of_node);
	while (1)
	{
		fscanf(fi, "%d%d", &p, &c);
		if (feof(fi)) break;

		if (max(p, c) >= G.size())
		{
			G.resize(max(p, c) + 1);
		}
		G[p].push_back(c);
	}
	fclose(fi);

	int nlink = 0;
	for (int i = 0; i < G.size(); i++)
	{
		nlink += G[i].size();
	}
	printf("LoadNetwork %s done, #node %d #link %d\n", fileName, number_of_node,  nlink);
}

void initilization(int number_of_node, vector<vector<int>> &vUser_vFriend_linkflag, vector<stStateOfInfection> &vUserState)
{
	//Initial infected state;
	stStateOfInfection ststate;

	vUserState.resize(number_of_node);
	for (int i = 0; i < number_of_node; i++)
	{
		ststate.userId = i;
		ststate.isInfected = false;//default susceptible
		ststate.infectedByWho = -1;//by nobody
		ststate.timestamp = -1;
		vUserState[i] = ststate;
	}
	printf("Initialize vUserState done\n");
	//Initial link flag
	vUser_vFriend_linkflag.resize(number_of_node);
	for (int i = 0; i < number_of_node; i++)
	{
		vUser_vFriend_linkflag[i].resize(number_of_node);
		for (int j = 0; j < number_of_node; j++)
		{
			vUser_vFriend_linkflag[i][j] = 0;
		}
	}
	printf("Initialize vUser_vFriend_linkflag done\n");
}

void record2rate(int T, float deltat)
{
	/*
	append the growth rate of each spreading instance to nodeAdd_rate and linkAdd_rate
	format: number per line for the growth rate @ time denoted by lineNumber.
	number
	*/
	FILE *fn, *fl;
	FILE *fno, *flo;

	int a, b;
	float t;
	map<float, int> day_node;
	map<float, int> day_link;
	map<float, int>::iterator it;

	int NT = ceil(T / deltat);
	for (int i = 1; i <= NT; i++)
	{
		day_node[i*deltat] = 0;
		day_link[i*deltat] = 0;
	}

	fn = fopen("nodeAdd_id_time.dat", "rb");
	while (1)
	{
		fscanf(fn, "%d%f", &a, &t);
		if (feof(fn)) break;
		day_node[t]++;
	}
	fno = fopen("nodeAdd_rate", "ab");
	for (it = day_node.begin(); it != day_node.end(); it++)
	{
		//fprintf(fno, "%d	%d\n", it->first, it->second);
		fprintf(fno, "%f %d	", it->first, it->second);
	}
	fprintf(fno, "\n");
	fclose(fno);
	fclose(fn);


	fl = fopen("linkAdd_ing_ed_time.dat", "rb");
	flo = fopen("linkAdd_rate", "ab");
	while (1)
	{
		fscanf(fl, "%d%d%f", &a, &b, &t);
		if (feof(fl)) break;
		day_link[t]++;
	}

	for (it = day_link.begin(); it != day_link.end(); it++)
	{
		//fprintf(flo, "%d	%d\n", it->first, it->second);
		fprintf(flo, "%f %d	", it->first, it->second);
	}
	fprintf(flo, "\n");
	fclose(fl);
	fclose(flo);
	printf("Dump generated rate file done.\n");

}


////PART2: Generators: NetTide-Process and NetTide-Survival
void NetTideProcess(vector<vector<int>> &vUser_vFriend, vector<vector<int>> &vUser_vFriend_linkflag, vector<stStateOfInfection> &vUserState, 
	double beta, double theta, double beta_prime, int t_max, float deltat, FILE *fo_node, FILE *fo_link)
{
	/*	
	NetTide-Node: node growth equation
		n-dot = beta / t^theta * n * (N - n )
	NetTide-Link: linking euqation
		e-dot = beta_prime / t^theta * n * ( alpha * (n-1)^gamma - e/n ) + 2*n-dot
	beta: growth rate of nodes
	theta: temporal fizzling exponent
	beta_prime: growth rate of links
	t_max: total time tick
	FILE *fo_node, *fo_link: dump evolving records.
		format:
			fo_node: node_id adding_time
			fo_link: parent_id child_id linking_time
	*/

	//1.Parameters given in main

	//2.Transforme grwoth rate to micro infection rate.
	double beta_micro, beta_prime_micro;
	int N = vUser_vFriend.size();    //number of nodes          
	float t = deltat;                        //time starts from minimum unit.

	double average_k = 0;             //avergae number of neighborhoods
	for (int i = 0; i < vUser_vFriend.size(); i++)
	{
		average_k += vUser_vFriend[i].size();
	}
	average_k /= vUser_vFriend.size();

	beta_micro = beta * N / average_k;
	beta_prime_micro = beta_prime;

	double rnd_number; //random number for probability infection
	vector<int> vInfectedUser;//operational target: infected user, for traverse

	//3. Set inital infected (existing) nodes 
	//default 0 th node is existing node in network. Infected.
	int seed_id;
	//seed_id = 0;
	seed_id = rand() % N;
	printf("Randomly choose infected seed, seed_id: %d\n", seed_id);
	stStateOfInfection ststate;
	ststate.userId = seed_id;
	ststate.infectedByWho = -1;
	ststate.isInfected = true;
	ststate.timestamp = t;

	vUserState[seed_id] = ststate;
	vInfectedUser.push_back(seed_id);
	fprintf(fo_node, "%d	%f\n", seed_id, t);
	t += deltat;

	int i, j;//for interation	
	queue<stStateOfInfection> queueInfectionForUpdata;
	double fizzling_kernel;
	//4. Begin NetTide process
	while (t <= t_max)
	{
		int infected_u, susceptible_u, infected_u_beLinked; //id of current infected user , susceptible, infected user to be Linked
		int n_friend, n_friendForLink, n_infected_now;  //number of friend of one user; number of infected user to be traversed

		n_infected_now = vInfectedUser.size();//in this tick, we only traver n_infected_now infected, push_back() have no effect.
		for (i = 0; i < n_infected_now; i++)
		{//For each infected user infected_u
			infected_u = vInfectedUser[i];

			//Step1: Try to infect infected_u's susceptible neighbors in G_0 with beta_micro.
			n_friend = vUser_vFriend[infected_u].size();
			for (j = 0; j < n_friend; j++)
			{
				susceptible_u = vUser_vFriend[infected_u][j];
				if (susceptible_u == infected_u)
				{//check self-loop, though we try to keep the self-loop not in the graph.
					//is self-loop
					continue;
				}
				if (vUserState[susceptible_u].isInfected == false)
				{
					//n-dot = beta / t^theta * (N - n ) * n 
					rnd_number = distribution(generator);
					fizzling_kernel = 1.0 / pow(t, theta);
					if (rnd_number <= beta_micro * fizzling_kernel * deltat)
					{
						ststate.userId = susceptible_u;
						ststate.infectedByWho = infected_u;
						ststate.isInfected = true;
						ststate.timestamp = t;
						// For updating all the infected state in Step3.
						queueInfectionForUpdata.push(ststate);


						//dump node adding and link adding records.
						fprintf(fo_node, "%d	%f\n", susceptible_u, t);

						//by infecting susceptible, there is a bidirection link at the same time.
						vUser_vFriend_linkflag[infected_u][susceptible_u] = 1;
						fprintf(fo_link, "%d	%d	%f\n", infected_u, susceptible_u, t);

						//comment out following two sentence to consider the reciprocity
						vUser_vFriend_linkflag[susceptible_u][infected_u] = 1;
						fprintf(fo_link, "%d	%d	%f\n", susceptible_u, infected_u, t);

					}
				}
			}//end step1

			//Step2: Try to link infected yet not connected neighbors in G_0 with beta_prime_micro
			n_friendForLink = vUser_vFriend[infected_u].size();
			for (j = 0; j < n_friendForLink; j++)
			{
				susceptible_u = vUser_vFriend[infected_u][j];
				if (susceptible_u == infected_u)
				{//check self-loop, though we try to keep the self-loop not in the graph.
					//is self-loop
					continue;
				}
				if (vUserState[susceptible_u].isInfected == true)
				{//try to build links between infected.
					//Link build between infected-infected and infected-newly infected:
					//J-dot = beta_prime * ( I * alpha * (I-1)^gamma - J ) + 2*I-dot
					infected_u_beLinked = susceptible_u;

					if (vUser_vFriend_linkflag[infected_u][infected_u_beLinked] != 1)
					{//Links between infected have two directions, 
						//like A-B = A->B + B->A, there are reciprocity.  
						rnd_number = distribution(generator);
						fizzling_kernel = 1.0 / pow(t, theta);
						if (rnd_number <= beta_prime_micro * fizzling_kernel * deltat)
						{
							vUser_vFriend_linkflag[infected_u][infected_u_beLinked] = 1;
							//dump node adding and link adding records.

							fprintf(fo_link, "%d	%d	%f\n", infected_u, infected_u_beLinked, t);

							//comment out following two sentence to consider the reciprocity
							vUser_vFriend_linkflag[infected_u_beLinked][infected_u] = 1;
							fprintf(fo_link, "%d	%d	%f\n", infected_u_beLinked, infected_u, t);


						}
					}
				}
			}//end step2

			// Step3: Update state of the infected for next time tick;
			while (!queueInfectionForUpdata.empty())
			{
				ststate = queueInfectionForUpdata.front();
				queueInfectionForUpdata.pop();
				vUserState[ststate.userId] = ststate;
				vInfectedUser.push_back(ststate.userId);
			}//end step3

		}//end traverse infected user
		t += deltat;
	}//end time tick.
	fflush(fo_node);
	fflush(fo_link);
}


void NetTideSurvival(int number_of_node, double beta, double theta, double beta_prime, double alpha, double gamma, int t_max, float deltat, FILE *fo_node, FILE *fo_link)
{
	/*
	NetTide-Node: node growth equation
	n-dot = beta / t^theta * n * (N - n )
	NetTide-Link: linking euqation
	e-dot = beta_prime / t^theta * n * ( alpha * (n-1)^gamma - e/n ) + 2*n-dot
	beta: growth rate of nodes
	theta: temporal fizzling exponent
	beta_prime: growth rate of links
	alpha: linear sparsity coefficient
	gamma: scaling sparsity exponent
	t_max: total time tick
	FILE *fo_node, *fo_link: dump evolving records.
	format:
	fo_node: node_id adding_time
	fo_link: parent_id child_id linking_time
	*/

	//1. Decalre. Parameters are initialized in main func.
	int N = number_of_node;    //number of nodes          
	float t = deltat;                        //time starts from minimum unit.
	double rnd_number; //random number for probability infection
	int i, j;//for interation	
	double fizzling_kernel;
	int nodesofar = 1, linksofar = 0;

	//2. Set inital infected (existing) nodes 
	//default 0 th node is existing node in network. Infected.
	int seed_id;
	seed_id = 0;
	printf("Randomly choose infected seed, seed_id: %d\n", seed_id);
	fprintf(fo_node, "%d	%f\n", seed_id, t);
	t += deltat;

	//3. Begin NetTide Survival process
	while (t <= t_max)
	{
		////int infected_u, susceptible_u, infected_u_beLinked; //id of current infected user , susceptible, infected user to be Linked
		int  n_infected_now, n_link_now, n_friendForLink;  //number of friend of one user; number of infected user to be traversed
		n_infected_now = nodesofar; //vInfectedUser.size();//in this tick, we only traver n_infected_now infected, push_back() have no effect.
		n_link_now = linksofar;                            //existing links so far.
		//n_friendForLink                                  //not established links to be linked. Links pool.

		//Step1: Try to infect  susceptible neighbors in survior pool. no network structure.
		for (i = 0; i < N - n_infected_now; i++)
		{//number of survivor pool: N - n_infected_now.
			rnd_number = distribution(generator);
			fizzling_kernel = 1.0 / pow(t, theta);
			if (rnd_number <= beta * fizzling_kernel * deltat * n_infected_now)
			{
				nodesofar++;
				//dump node adding and link adding records.
				fprintf(fo_node, "%d	%f\n", nodesofar, t);
				//by infecting susceptible, there is a bidirection link at the same time.
				linksofar++;
				fprintf(fo_link, "%d	%d	%f\n", -1, nodesofar, t);			
				//comment out following two sentence to consider the reciprocity
				linksofar++;
				fprintf(fo_link, "%d	%d	%f\n", nodesofar, -1, t);
			}
		}//End step1

		//Step2: Try to link infected yet not connected neighbors 
		n_friendForLink = alpha *  n_infected_now *   pow((n_infected_now - 1), gamma) - n_link_now;//vUser_vFriend[infected_u].size();
		for (j = 0; j < n_friendForLink; j++)
		{
			rnd_number = distribution(generator);
			fizzling_kernel = 1.0 / pow(t, theta);
			if (rnd_number <= beta_prime * fizzling_kernel * deltat)
			{
				//vUser_vFriend_linkflag[infected_u][infected_u_beLinked] = 1;
				//dump node adding and link adding records.
				linksofar++;
				fprintf(fo_link, "%d	%d	%f\n", -1, linksofar, t);

				//comment out following two sentence to consider the reciprocity
				//vUser_vFriend_linkflag[infected_u_beLinked][infected_u] = 1;
				linksofar++;
				fprintf(fo_link, "%d	%d	%f\n", -1, linksofar, t);
			}
		}//end step2

		t += deltat;
	}//end time tick.

	fflush(fo_node);
	fflush(fo_link);
}


////PART3: How to use Generators.
void run_NetTide_Process()
{
	/*
	Run NetTide Process in Section 3.4. Results in Section 4.4
	NetTide-Process captures the micro-level stochastic interactions within a network.
	*/

	time_t ltstart, ltend;
	ltstart = time(NULL);
	// Declare variables
	vector<vector<int>> vUser_vFriend; //(NUMBER_OF_NODE); underlying network structure for node adoption
	vector<vector<int>> vUser_vFriend_linkflag; //[NUMBER_OF_NODE][NUMBER_OF_NODE]; mark link: addjacency matrix of link;
	vector<stStateOfInfection> vUserState; //[NUMBER_OF_NODE]; mark user state, like hash_table
	FILE *fo_node;
	FILE *fo_link;

	//Step 1, Set the total number of nodes.
	int number_of_node = 10000;
	
	//Step 2, Load the undelrying network (organizational structure).
	loadNetwork(number_of_node, vUser_vFriend, "complete_10000_0.01_V2.txt");

	/* The output files, documenting the rate over time.
	   each row corresponds to one instance, e.g. 1st row of nodeAdd_rate and linkAdd_rate are from 1st instance 
	   Format: 
	          nodeAdd_rate: t1 n(t1) t2 n(t2) ...
	          linkAdd_rate: t1 e(t1) t2 e(t2) ...
	*/
	fopen("nodeAdd_rate", "wb"); //clear the file
	fopen("linkAdd_rate", "wb"); // clear the file
	
	//Step 3, Set the parameters for simulation, definitions in function NetTideProcess    
	int iterN = 10;  1001; 200; 1001;//the number of simulation times. 
	float deltat = 1;
	double beta = 2.485e-4;// 2.47e-4;//2.45e-4;;//2.5e-4;//1.28e-4;//2e-4;//0.5*1.28e-4; 0.0002;//0.00003; //0.0002
	double theta = 1;
	double beta_prime = 0.48;// 0.49;// 0.5; //0.53; //0.55;  //0.6; //0.7;
	int t_max = 1000; 500; 1000;//% maximum time-tick


	for (int i = 1; i <= iterN; i++)
	{
		printf("\n>>>>>BEGIN ITER: %d\n", i);
		fo_node = fopen("nodeAdd_id_time.dat", "wb");
		fo_link = fopen("linkAdd_ing_ed_time.dat", "wb");

		initilization(number_of_node, vUser_vFriend_linkflag, vUserState);

		NetTideProcess(vUser_vFriend, vUser_vFriend_linkflag, vUserState, beta, theta, beta_prime, t_max, deltat, fo_node, fo_link);
		//Only the last growth records are stored because of overwritting.
		//All the dynamics records are stored in nodeAdd_rate and linkAdd_rate.
		fclose(fo_node);
		fclose(fo_link);
		record2rate(t_max, deltat);
		ltend = time(NULL);
		printf("DONE ITER %d, cost time diff: %f\n", i, difftime(ltend, ltstart));
	}

}

void run_NetTide_Survival()
{
	/*
	Run NetTide Survival process in Section 3.4. Results in Section 4.4
	NetTide-Survival captures the growth dynamics from hazard rates (hazard process).
	*/
	time_t ltstart, ltend;
	ltstart = time(NULL);
	// Declare variables

	FILE *fo_node;
	FILE *fo_link;

	/* The output files, documenting the rate over time.
	   each row corresponds to one instance, e.g. 1st row of nodeAdd_rate and linkAdd_rate are from 1st instance 
	   Format: 
	          nodeAdd_rate: t1 n(t1) t2 n(t2) ...
	          linkAdd_rate: t1 e(t1) t2 e(t2) ...
	*/
	fopen("nodeAdd_rate", "wb"); //clear the file
	fopen("linkAdd_rate", "wb"); // clear the file

	//Step 1, Set the total number of nodes.
	int number_of_node = 10000;
	//Step 2, Set the parameters for simulation, definitions in function NetTideProcess   
	int iterN = 10;  1001; 200; 1001;//the number of simulation times. 
	float deltat = 1;
	double beta = 2.45e-4; 
	double theta = 1;
	double beta_prime = 0.5;// 0.49;// 0.5; //0.53; //0.55;  //0.6; //0.7;
	double alpha = 0.9;
	double gamma = 0.52;
	int t_max = 1000;  //% maximum time-tick


	for (int i = 1; i <= iterN; i++)
	{
		printf("\n>>>>>BEGIN ITER: %d\n", i);
		fo_node = fopen("nodeAdd_id_time.dat", "wb");
		fo_link = fopen("linkAdd_ing_ed_time.dat", "wb");

		NetTideSurvival(number_of_node, beta, theta, beta_prime, alpha, gamma, t_max, deltat, fo_node, fo_link);
		//Only the last growth records are stored because of overwritting.
		//All the dynamics records are stored in nodeAdd_rate and linkAdd_rate.
		fclose(fo_node);
		fclose(fo_link);
		record2rate(t_max, deltat);
		ltend = time(NULL);
		printf("DONE ITER %d, cost time diff: %f\n", i, difftime(ltend, ltstart));
	}

}

int main()
{
//Choose one generator to run.

	run_NetTide_Process();
	//run_NetTide_Survival();
}

