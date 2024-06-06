#include <iostream> 
#include <math.h> 
#include <fstream> 
#include <sstream> 
#include <stdlib.h>  
#include <random> 
#include <iomanip> 
#include "omp.h" 
#include <numeric> 

using namespace std; 

#define N pow(10,5) 
#define dt 1.0*pow(10,-4) 
#define sigma 5.0
#define s 2.0 
#define lambda_c 1.0
#define lambda_h 5.0  
#define kB 1.380649*pow(10,-5)
//#define kB 1.0
#define Th 363.15
//#define Th 1.0
#define Tc 278.15
//#define Tc 1.0

default_random_engine generator;
normal_distribution<double> distribution(0.0, sqrt(1.0*dt));

default_random_engine generatorY;
normal_distribution<double> distributionY(0.0, sqrt(pow(sigma,2)*s/2));

ofstream out6("./Brownian_x2_Y2_final_StirlingCycle_t0_tau0p01.dat");
ofstream out7("./Brownian_x2_Y2_final_StirlingCycle_ttauby2_tau0p01.dat");
ofstream out8("./Brownian_x2_Y2_final_StirlingCycle_ttau_tau0p01.dat");
ofstream out9("./Brownian_x2_Y2_finalcycle_StirlingCycle_tau0p01.dat");

vector<double> diffusion(vector<int> &index_Y, vector<double> Y, int size)
{
	vector<double> dummy_Y(size, 0);
	dummy_Y.reserve(index_Y.size());

	vector<double> deltaW(size, 0);
	vector<double> k1(size, 0);

	//double std_gauss = sigma*sqrt(dt);

	double noise_eta;	
	for(int i=0;i<size;i++)
	{	
		dummy_Y[i] = Y[i];
		noise_eta = distribution(generator);
		deltaW[i] = noise_eta;
		k1[i] = -dummy_Y[i]/s;
		dummy_Y[i] = dummy_Y[i]+dt*k1[i]+(sigma*deltaW[i]);
	}

	//out<<fixed<<setprecision(10)<<"   "<<noise_eta<<endl;
	return dummy_Y;
}

int main()
{
	int N_tau=1;
	double Avg_Q[N_tau][24][2];
	
		int n_tau = 0;
		double stept_cycle; 
				
		stept_cycle = 0+100*(n_tau+1);
		
		double t_cycle = stept_cycle*dt;
		
		vector<double> dof_for(N, 0);
		vector<double> dof_back(N, 0);
		vector<double> deltaW(N, 0);
		vector<double> dof_for_tauby2(N, 0);
		vector<double> dof_initial(N, 0);
		vector<double> cum_dof(N, 0);

		double std_x2 = 1, mean_x2 = 0, mean_x2_prev = 0, mean_x2_cum = 0, precision1 = 0.000001, precision = 0.0007;

		cout<<"start = "<<0<<", t_cycle = "<<t_cycle<<", stept_cycle = "<<stept_cycle<<endl;
		
		int n_cycle = 0;
		int CLM_start = 0;
		int enter = 0;
		int n_cycle_TPSS = 0;
		int n_cycle_stats = 0;
        int times=1;

		vector<double> Y(N, 0);
		vector<int> index_Y(N);

		for(int i=0; i<N; i++)
		{
			Y[i] = distributionY(generatorY);		
		}		

		while(enter != 2)
		{
			vector<double> w1(N, 0);
			vector<double> w3(N, 0);
			vector<double> x_back(stept_cycle/2, 0);
			vector<double> x2_back(stept_cycle/2, 0);
			vector<double> x3_back(stept_cycle/2, 0);
			vector<double> x4_back(stept_cycle/2, 0);
			vector<double> Y_back(stept_cycle/2, 0);
			vector<double> Y2_back(stept_cycle/2, 0);

			mean_x2 = 0;

			for(int stept=0;stept<=stept_cycle/2;stept++)
			{		    
				double t = stept*dt;
				double noise_xi;	
				
				//if(stept%1000==0)
				//{cout<<"stept = "<<stept<<endl;}
				
				if(n_cycle>0 || stept!=0)	
				{Y = diffusion(index_Y, Y, N);}

				double x_mean=0, x2_mean=0, x3_mean=0, x4_mean=0, Y_mean=0, Y2_mean=0;
				for(int i=0;i<N;i++)
				{	
					int step = 0;
					noise_xi = distribution(generator);
					deltaW[i] = noise_xi;
					if(stept==0)
					{
						//dof_for[i] = distributionx(generatorx);
						dof_initial[i] = dof_for[i];
					}

					if(enter==1)
					{
						if(stept==0)
						{
							out6<<fixed<<setprecision(10)<<stept<<"  "<<dof_for[i]<<"  "<<Y[i]<<endl;
						}
						if(stept==stept_cycle/2)
						{
							out7<<fixed<<setprecision(10)<<stept<<"  "<<dof_for[i]<<"  "<<Y[i]<<endl;
						}

						x_mean = x_mean+dof_for[i];
						x2_mean = x2_mean+pow(dof_for[i],2);
						x3_mean = x3_mean+pow(dof_for[i],3);
						x4_mean = x4_mean+pow(dof_for[i],4);
						Y_mean = Y_mean+Y[i];
						Y2_mean = Y2_mean+pow(Y[i],2);

						if(i==N-1)
						{out9<<fixed<<setprecision(10)<<t<<"  "<<(lambda_h-(2.0*(lambda_h-lambda_c)*1.0*t/(1.0*t_cycle)))<<"  "<<1.0*x_mean/N<<"  "<<1.0*x2_mean/N<<"  "<<1.0*x3_mean/N<<"  "<<1.0*x4_mean/N<<"  "<<1.0*Y_mean/N<<"  "<<1.0*Y2_mean/N<<endl;}
					}

					dof_for[i] = dof_for[i]-((lambda_h-(2.0*(lambda_h-lambda_c)*1.0*t/(1.0*t_cycle)))*dof_for[i]*pow(1.0*Y[i],2)*dt)+(sqrt(2.0*kB*Th*pow(1.0*Y[i],2))*deltaW[i]);
					//dof_for[i] = dof_for[i]-((lambda_c-((1.0-2.0*(1.0*t/1.0*t_cycle))*(lambda_h-lambda_c)))*dof_for[i]*pow(1.0*Y[i],2)*dt)+(sqrt(2.0*kB*Tc*pow(1.0*Y[i],2))*deltaW[i]);
					
					// if(i==N-1 && (stept >40))
					// {
					// 	// cout<<"At Th = "<<endl;
					// 	// cout<<dof_for[i]<<"   "<<pow(1.0*Y[i],2)<<"   "<<sqrt(2.0*kB*Th*pow(1.0*Y[i],2))*deltaW[i]<<"   "<<((lambda_h-(2.0*(lambda_h/lambda_c)*1.0*t/1.0*t_cycle))*dof_for[i]*pow(1.0*Y[i],2)*dt)<<endl;
					// 	cout<<"At Tc = "<<endl;
					// 	cout<<dof_for[i]<<"   "<<pow(1.0*Y[i],2)<<"   "<<sqrt(2.0*kB*Tc*pow(1.0*Y[i],2))*deltaW[i]<<"   "<<((lambda_c-((1.0-2.0*(t/t_cycle))*(lambda_h-lambda_c)))*dof_for[i]*pow(1.0*Y[i],2)*dt)<<endl;
					// }	

					w1[i] = w1[i] + pow(dof_for[i],2);

					//cout<<t<<"   "<<1.0*t_cycle/2<<endl;
					if(stept==(stept_cycle/2))
					{
						int t_stept = (stept_cycle/2)+1;
						double t_step = dt*((stept_cycle/2)+1);
						dof_back[i] = dof_for[i];

						vector<double> Y1(1, 0);
						vector<int> index_Y1(1);
						Y1[0] = Y[i];

						while(t_stept<=stept_cycle)
						{
							// if(i==0)
							// {
								// if(t_stept%1000==0)
								// {cout<<"t_stept = "<<t_stept<<endl;}
								Y1 = diffusion(index_Y1, Y1, 1);
							//}
					
							noise_xi = distribution(generator);
							deltaW[i] = noise_xi;
							dof_back[i] = dof_back[i]-((lambda_c-((1.0-2.0*(1.0*t_step/(1.0*t_cycle)))*(lambda_h-lambda_c)))*dof_back[i]*pow(1.0*Y1[0],2)*dt)+(sqrt(2.0*kB*Tc*pow(1.0*Y1[0],2))*deltaW[i]);
							// if(i==N-1 && (t_stept > 90))
							// {
							// 	cout<<"At Tc = "<<endl;
							// 	cout<<dof_back[i]<<"   "<<pow(1.0*Y[i],2)<<"   "<<sqrt(2.0*kB*Tc*pow(1.0*Y[i],2))*deltaW[i]<<"   "<<((lambda_c-((1.0-2.0*(t_step/t_cycle))*(lambda_h-lambda_c)))*dof_back[i]*pow(1.0*Y[i],2)*dt)<<endl;
							// }	

							w3[i] = w3[i] + pow(dof_back[i],2);	
							
							if(t_stept==stept_cycle)
							{
								dof_for_tauby2[i] = dof_for[i];
								dof_for[i] = dof_back[i];
								mean_x2 = mean_x2 + pow(dof_back[i],2);
								if(CLM_start==1)
								{
								    mean_x2_cum = mean_x2_cum + pow(dof_back[i],2);
								    cum_dof[i] = cum_dof[i] + pow(dof_back[i],2);
                                }   
							}

							if(enter==1)
							{
								if(t_stept==stept_cycle)
								{
									out8<<fixed<<setprecision(10)<<t_stept<<"  "<<dof_back[i]<<"  "<<Y[i]<<endl;
								}
								x_back[step] = x_back[step] + dof_back[i];
								x2_back[step] = x2_back[step] + pow(dof_back[i],2);
								x3_back[step] = x3_back[step] + pow(dof_back[i],3);
								x4_back[step] = x4_back[step] + pow(dof_back[i],4);
								Y_back[step] = Y_back[step] + Y1[0];
								Y2_back[step] = Y2_back[step] + pow(Y1[0],2);

								if(i==N-1)
								{
									out9<<fixed<<setprecision(10)<<t_step<<"  "<<(lambda_c-((1.0-2.0*(1.0*t_step/(1.0*t_cycle)))*(lambda_h-lambda_c)))<<"  "<<1.0*x_back[step]/N<<"  "<<1.0*x2_back[step]/N<<"  "<<1.0*x3_back[step]/N<<"  "<<1.0*x4_back[step]/N<<"  "<<1.0*Y_back[step]/N<<"  "<<1.0*Y2_back[step]/N<<endl;	
								}

								step++;
							}

							t_stept++;
							t_step = dt*t_stept;
						}
						
						Y[i] = Y1[0];
					}
				}
			}	

			if(CLM_start==1)
			{
				std_x2=0;
				for(int i=0;i<N;i++)
				{
					std_x2 = std_x2 + pow(((cum_dof[i]/times)-(mean_x2_cum/(N*times))),2);
				}
				std_x2 = sqrt(std_x2/N);
                // if(enter==0)
				// {
				// 	cout<<"n_cycle_TPSS = "<<n_cycle_TPSS<<", std_x2 = "<<std_x2<<", times = "<<times<<", n_cycle = "<<n_cycle<<", 1/tau = "<<int(1.0/t_cycle)<<endl;	
				// }	
				// else
				// {
				// 	cout<<"n_cycle_stats = "<<n_cycle_stats<<", std_x2 = "<<std_x2<<", times = "<<times<<", n_cycle = "<<n_cycle<<", 1/tau = "<<int(1.0/t_cycle)<<endl;	
				// }
				times++;
			}

            n_cycle++;	

			if(1.0*fabs(mean_x2-mean_x2_prev)/N>precision1 && CLM_start==0)
			{
				n_cycle_TPSS++;
				cout<<"n_cycle_TPSS = "<<n_cycle_TPSS<<", std_x2 = "<<std_x2<<endl;	
				cout<<"n_cycle_TPSS = "<<n_cycle_TPSS<<", deviation = "<<1.0*fabs(mean_x2-mean_x2_prev)/N<<endl;	
                mean_x2_prev = mean_x2;
			}	
			else
			{
				if(CLM_start==0)
				{cout<<"n_cycle_TPSS = "<<n_cycle_TPSS<<", deviation = "<<1.0*fabs(mean_x2-mean_x2_prev)/N<<endl;}
				n_cycle_stats++;
				cout<<"n_cycle_stats = "<<n_cycle_stats<<", std_x2 = "<<std_x2<<endl;	
				CLM_start = 1;	
			}

			if(CLM_start==1)
			{
				
				if(std_x2<precision)
				{		
					enter++;
				}

				//cout<<"n_cycle_stats = "<<"   "<<Cyc_for_TPSS-n_cycle_TPSS+1<<endl;
			}

			cout<<"t_cycle = "<<t_cycle<<", n_cycle = "<<n_cycle<<endl;	
		}	
}
