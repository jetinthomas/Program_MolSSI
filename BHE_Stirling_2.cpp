#include <iostream> 
#include <math.h> 
#include <fstream> 
#include <sstream> 
#include <stdlib.h>  
#include <random> 
#include <iomanip> 
#include <cmath>
#include "omp.h" 
#include <numeric> 

using namespace std; 

#define N pow(10,5) 
#define dt 1.0*pow(10,-4) 
#define sigma 1.0
#define tau 1.0 
#define eta1 0.004
#define a 0.5
#define gma 6.0*M_PI*eta1*a
//#define gma 2
#define lambda_c 1.0
#define lambda_h 5.0  
#define kB 1.380649*pow(10,-5)
#define Th 363.15
#define Tc 278.15

default_random_engine generator;
normal_distribution<double> distribution(0.0, sqrt(1.0*dt));

//default_random_engine generatorx;
//normal_distribution<double> distributionx(0.0, sqrt(1.0*kB*Th/lambda_h));

ofstream out("./Brownian_Stirling_AvgQnt_with_time_check.dat");
ofstream out1("./Brownian_Stirling_tau0p01.dat");
// ofstream out2("./Brownian_Stirling_tau7.dat");
// ofstream out3("./Brownian_Stirling_tau70.dat");
ofstream out4("./Brownian_Dispsq_with_StirlingCycles_check.dat");
ofstream out5("./Brownian_Stirling_Error_with_time_check.dat");

int main()
{
	int N_tau=1;
	double Avg_Q[N_tau][24][2];
	//omp_set_num_threads(N_tau);

	//cout<<gma<<"   "<<kB<<"   "<<sqrt((2.0*kB*Th)/(1.0*gma))<<"   "<<sqrt(1.0*kB*Th/lambda_h)<<"   "<<lambda_h*dt/gma<<endl;	

	
		int n_tau = 0;
		double stept_cycle = 0 + 1000*(n_tau+1);
		double t_cycle = stept_cycle*dt;
		
		//vector<double> Y(N, 1);
		//vector<int> index_Y(N);
		vector<double> dof_for(N, 0);
		vector<double> dof_back(N, 0);
		vector<double> deltaW(N, 0);
		vector<double> dof_for_tauby2(N, 0);
		vector<double> dof_initial(N, 0);
		vector<double> cum_dof(N, 0);

		//for(int i=0;i<N;i++)
		//{
		//	dof_for[i] = distributionx(generatorx);
		//}
	
		double std_x2 = 1, mean_x2 = 0, mean_x2_prev = 0, mean_x2_cum = 0, precision1 = 0.000001, precision = 0.0007;

		cout<<"start = "<<0<<", t_cycle = "<<t_cycle<<", stept_cycle = "<<stept_cycle<<endl;
			
		int n_cycle = 0;
		int CLM_start = 0;
		int enter = 0;
		int n_cycle_TPSS = 0;
		int n_cycle_stats = 0;
        int times=1;

		while(enter != 2)
		{
			vector<double> w1(N, 0);
			vector<double> w3(N, 0);
			vector<double> x2_TPSS_back(stept_cycle/2, 0);
			vector<double> mean_TPSS_back(stept_cycle/2, 0);
            mean_x2 = 0;

			for(int stept=0;stept<=stept_cycle/2;stept++)
			{		    
				double t = stept*dt;
				//Y = diffusion(index_Y, Y);
				double noise_xi;	
				double x2_TPSS = 0;
				double mean_TPSS = 0;
				int ind = 0;
								

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

					dof_for[i] = dof_for[i]-((lambda_h-(2*(lambda_h-lambda_c)*t/t_cycle))*dof_for[i]*dt/(1.0*gma))+(sqrt(2.0*kB*Th/(1.0*gma))*deltaW[i]);

					//if(i==0 && (stept%10)==0)
					//{cout<<"t = "<<t<<", dof_for = "<<pow(dof_for[i],2)<<endl;}
				
					w1[i] = w1[i] + pow(dof_for[i],2);

					//cout<<t<<"   "<<1.0*t_cycle/2<<endl;
					if(stept==stept_cycle/2)
					{
						int t_stept = (stept_cycle/2)+1;
						double t_step = dt*((stept_cycle/2)+1);
						dof_back[i] = dof_for[i];
						while(t_stept<=stept_cycle)
						{
							noise_xi = distribution(generator);
							deltaW[i] = noise_xi;
							dof_back[i] = dof_back[i]-((lambda_c-((1-(2*(t_step/t_cycle)))*(lambda_h-lambda_c)))*dof_back[i]*dt/(1.0*gma))+(sqrt(2.0*kB*Tc/(1.0*gma))*deltaW[i]);
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

							//if(i==0 && (t_stept%10)==0)
							//{
							//	cout<<"i = "<<i<<", dof_tauby2 = "<<pow(dof_for_tauby2[i],2)<<", dof_back = "<<pow(dof_back[i],2)<<endl;							
							//}
						
							if((t_stept+1)%10==0)
							{
								x2_TPSS_back[step] = x2_TPSS_back[step] + pow(dof_back[i],2);
								mean_TPSS_back[step] = mean_TPSS_back[step] + dof_back[i];
								//cout<<"t_step = "<<t_stept<<", t_step = "<<t_step<<", step = "<<step<<", dof_back = "<<pow(dof_back[i],2)<<endl;
								if(i==N-1)
								{
								 	out4<<fixed<<setprecision(10)<<"   "<<t_cycle<<"   "<<n_cycle<<"   "<<t_step+dt<<"   "<<sqrt((1.0*x2_TPSS_back[step]/N)-pow((mean_TPSS_back[step]/N),2))<<endl;
								}
								step++;
							}

							t_stept++;
							t_step = dt*t_stept;
						}
					}
						
					if((stept+1)%10==0 && stept!=(stept_cycle/2))
					{
						ind = 1.0;
						x2_TPSS = x2_TPSS + pow(dof_for[i],2);
						mean_TPSS = mean_TPSS + dof_for[i];
						//cout<<"stept = "<<stept<<", x2_TPSS = "<<x2_TPSS<<", dof_for = "<<dof_for[i]<<endl;
					}
				}

				if(ind==1)
				{out4<<fixed<<setprecision(10)<<"   "<<t_cycle<<"   "<<n_cycle<<"   "<<t+dt<<"   "<<sqrt((1.0*x2_TPSS/N)-pow(mean_TPSS/N,2))<<endl;}
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

				double mean_w=0, mean_w2=0, mean_w3=0, mean_w4=0, 
				mean_q1=0, mean_q1_2=0, mean_q1_3=0, mean_q1_4=0,
				mean_q2=0, mean_q2_2=0, mean_q2_3=0, mean_q2_4=0,
				mean_dU=0, mean_P=0, mean_P2=0, mean_P3=0, mean_P4=0, 
                mean_eta=0, mean_eta2=0, mean_eta3=0, mean_eta4=0;

				double eta_bar;
				
				for(int i=0;i<N;i++)
				{	   
					double du1, q1, du2, w2, du3, q2, du4, w4, w, dU, P, eta;
					
					du1 = ((1.0*lambda_c/2)*pow(dof_for_tauby2[i],2)) - ((1.0*lambda_h/2)*pow(dof_initial[i],2));
					q1 = ((-1.0*(lambda_h-lambda_c)*dt/t_cycle)*w1[i]) - du1;
				    du2 = 0;
				   	w2 = 0;
				    du3 = ((1.0*lambda_h/2)*pow(dof_back[i],2)) - ((1.0*lambda_c/2)*pow(dof_for_tauby2[i],2));
				    q2 = ((1.0*(lambda_h-lambda_c)*dt/t_cycle)*w3[i]) - du3;
				    du4 = 0;
				    w4 = 0;

				    w = ((-1.0*(lambda_h-lambda_c)*dt/t_cycle)*w1[i])+w2+((1.0*(lambda_h-lambda_c)*dt/t_cycle)*w3[i])+w4;
				    dU = du1+du2+du3+du4;
				    P = -1.0*w/t_cycle;
				    eta = 1.0*w/q1;	

					//cout<<w<<"   "<<q1<<"   "<<eta<<endl;

					mean_w = mean_w + w;
					mean_w2 = mean_w2 + pow(w,2);
					mean_w3 = mean_w3 + pow(w,3);
					mean_w4 = mean_w4 + pow(w,4);

					mean_q1 = mean_q1 + q1;
					mean_q1_2 = mean_q1_2 + pow(q1,2);
					mean_q1_3 = mean_q1_3 + pow(q1,3);
					mean_q1_4 = mean_q1_4 + pow(q1,4);

					mean_q2 = mean_q2 + q2;
					mean_q2_2 = mean_q2_2 + pow(q2,2);
					mean_q2_3 = mean_q2_3 + pow(q2,3);
					mean_q2_4 = mean_q2_4 + pow(q2,4);

					mean_dU = mean_dU + dU;

					mean_P = mean_P + P;
                    mean_P2 = mean_P2 + pow(P,2);
                    mean_P3 = mean_P3 + pow(P,3);
                    mean_P4 = mean_P4 + pow(P,4);

					
					mean_eta = mean_eta + eta;
                    mean_eta2 = mean_eta2 + pow(eta,2);
                    mean_eta3 = mean_eta3 + pow(eta,3);
                    mean_eta4 = mean_eta4 + pow(eta,4);
						
					if(enter==2)
					{out1<<fixed<<setprecision(10)<<stept_cycle*dt<<"   "<<w<<"  "<<q1<<"  "<<q2<<"  "<<dU<<"  "<<P<<"   "<<eta<<"   "<<dof_back[i]<<endl;}
	
					//if(stept_cycle==1000)
					//{
					//	out1<<fixed<<setprecision(10)<<stept_cycle*dt<<"   "<<w<<"   "<<eta<<endl;
					//}
					// if(stept_cycle==35000)
					// {
					// 	out2<<fixed<<setprecision(10)<<stept_cycle*dt<<"   "<<w<<"   "<<eta<<endl;		
					// }
					// if(stept_cycle==350000)
					// {
					// 	out3<<fixed<<setprecision(10)<<stept_cycle*dt<<"   "<<w<<"   "<<eta<<endl;
					// }
				}	

				mean_w = mean_w/N;
				mean_w2 = mean_w2/N; 
				mean_w3 = mean_w3/N;
				mean_w4 = mean_w4/N;

			    mean_q1 = mean_q1/N;
				mean_q1_2 = mean_q1_2/N;
				mean_q1_3 = mean_q1_3/N;
				mean_q1_4 = mean_q1_4/N;

			    mean_q2 = mean_q2/N;
				mean_q2_2 = mean_q2_2/N;
				mean_q2_3 = mean_q2_3/N;
				mean_q2_4 = mean_q2_4/N;

			    mean_dU = mean_dU/N;

			    mean_P = mean_P/N;
                mean_P2 = mean_P2/N;
                mean_P3 = mean_P3/N;
                mean_P4 = mean_P4/N;

			    mean_eta = mean_eta/N;
                mean_eta2 = mean_eta2/N;
                mean_eta3 = mean_eta3/N;
                mean_eta4 = mean_eta4/N;

			    eta_bar = 1.0*mean_w/mean_q1;

                Avg_Q[n_tau][0][0] = t_cycle;
				Avg_Q[n_tau][0][1] = t_cycle;

			    Avg_Q[n_tau][1][0] = Avg_Q[n_tau][1][0] + mean_w;
				Avg_Q[n_tau][1][1] = Avg_Q[n_tau][1][1] + pow(mean_w,2);
				Avg_Q[n_tau][2][0] = Avg_Q[n_tau][2][0] + mean_w2;
				Avg_Q[n_tau][2][1] = Avg_Q[n_tau][2][1] + pow(mean_w2,2);
				Avg_Q[n_tau][3][0] = Avg_Q[n_tau][3][0] + mean_w3;
				Avg_Q[n_tau][3][1] = Avg_Q[n_tau][3][1] + pow(mean_w3,2);
				Avg_Q[n_tau][4][0] = Avg_Q[n_tau][4][0] + mean_w4;
				Avg_Q[n_tau][4][1] = Avg_Q[n_tau][4][1] + pow(mean_w4,2);

			    Avg_Q[n_tau][5][0] = Avg_Q[n_tau][5][0] + mean_q1;
				Avg_Q[n_tau][5][1] = Avg_Q[n_tau][5][1] + pow(mean_q1,2);
				Avg_Q[n_tau][6][0] = Avg_Q[n_tau][6][0] + mean_q1_2;
				Avg_Q[n_tau][6][1] = Avg_Q[n_tau][6][1] + pow(mean_q1_2,2);
				Avg_Q[n_tau][7][0] = Avg_Q[n_tau][7][0] + mean_q1_3;
				Avg_Q[n_tau][7][1] = Avg_Q[n_tau][7][1] + pow(mean_q1_3,2);
				Avg_Q[n_tau][8][0] = Avg_Q[n_tau][8][0] + mean_q1_4;
				Avg_Q[n_tau][8][1] = Avg_Q[n_tau][8][1] + pow(mean_q1_4,2);

			    Avg_Q[n_tau][9][0] = Avg_Q[n_tau][9][0] + mean_q2;
				Avg_Q[n_tau][9][1] = Avg_Q[n_tau][9][1] + pow(mean_q2,2);
				Avg_Q[n_tau][10][0] = Avg_Q[n_tau][10][0] + mean_q2_2;
				Avg_Q[n_tau][10][1] = Avg_Q[n_tau][10][1] + pow(mean_q2_2,2);
				Avg_Q[n_tau][11][0] = Avg_Q[n_tau][11][0] + mean_q2_3;
				Avg_Q[n_tau][11][1] = Avg_Q[n_tau][11][1] + pow(mean_q2_3,2);
				Avg_Q[n_tau][12][0] = Avg_Q[n_tau][12][0] + mean_q2_4;
				Avg_Q[n_tau][12][1] = Avg_Q[n_tau][12][1] + pow(mean_q2_4,2);

			    Avg_Q[n_tau][13][0] = Avg_Q[n_tau][13][0] + mean_dU;
				Avg_Q[n_tau][13][1] = Avg_Q[n_tau][13][1] + pow(mean_dU,2);
					
			    Avg_Q[n_tau][14][0] = Avg_Q[n_tau][14][0] + mean_P;
				Avg_Q[n_tau][14][1] = Avg_Q[n_tau][14][1] + pow(mean_P,2);
                Avg_Q[n_tau][15][0] = Avg_Q[n_tau][15][0] + mean_P2;
				Avg_Q[n_tau][15][1] = Avg_Q[n_tau][15][1] + pow(mean_P2,2);
                Avg_Q[n_tau][16][0] = Avg_Q[n_tau][16][0] + mean_P3;
				Avg_Q[n_tau][16][1] = Avg_Q[n_tau][16][1] + pow(mean_P3,2);
                Avg_Q[n_tau][17][0] = Avg_Q[n_tau][17][0] + mean_P4;
				Avg_Q[n_tau][17][1] = Avg_Q[n_tau][17][1] + pow(mean_P4,2);
					
					
				Avg_Q[n_tau][18][0] = Avg_Q[n_tau][18][0] + mean_eta;
				Avg_Q[n_tau][18][1] = Avg_Q[n_tau][18][1] + pow(mean_eta,2);
                Avg_Q[n_tau][19][0] = Avg_Q[n_tau][19][0] + mean_eta2;
				Avg_Q[n_tau][19][1] = Avg_Q[n_tau][19][1] + pow(mean_eta2,2);
                Avg_Q[n_tau][20][0] = Avg_Q[n_tau][20][0] + mean_eta3;
				Avg_Q[n_tau][20][1] = Avg_Q[n_tau][20][1] + pow(mean_eta3,2);
                Avg_Q[n_tau][21][0] = Avg_Q[n_tau][21][0] + mean_eta4;

				Avg_Q[n_tau][21][1] = Avg_Q[n_tau][21][1] + pow(mean_eta4,2);

			    Avg_Q[n_tau][22][0] = Avg_Q[n_tau][22][0] + eta_bar;
				Avg_Q[n_tau][22][1] = Avg_Q[n_tau][22][1] + pow(eta_bar,2);

				Avg_Q[n_tau][23][0] = n_cycle_stats;
				Avg_Q[n_tau][23][1] = n_cycle_stats;
				

				if(eta_bar>1)
				{
					cout<<"error = "<<n_tau<<"   "<<eta_bar<<"   "<<mean_w<<"   "<<mean_q1<<endl;
				}

				//cout<<"n_cycle_stats = "<<"   "<<Cyc_for_TPSS-n_cycle_TPSS+1<<endl;
			}

			cout<<"t_cycle = "<<t_cycle<<", n_cycle = "<<n_cycle<<endl;	
		}
			
	
	for(int n_tau=0; n_tau<N_tau; n_tau++)
	{
		int k=0;
		out<<fixed<<setprecision(10)<<Avg_Q[n_tau][0][k]<<"   "<<(Avg_Q[n_tau][1][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][2][k]/Avg_Q[n_tau][23][k])
        	<<"   "<<(Avg_Q[n_tau][3][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][4][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][5][k]/Avg_Q[n_tau][23][k])
		 <<"   "<<(Avg_Q[n_tau][6][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][7][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][8][k]/Avg_Q[n_tau][23][k])
		 <<"   "<<(Avg_Q[n_tau][9][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][10][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][11][k]/Avg_Q[n_tau][23][k])
		 <<"   "<<(Avg_Q[n_tau][12][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][13][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][14][k]/Avg_Q[n_tau][23][k])
		 <<"   "<<(Avg_Q[n_tau][15][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][16][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][17][k]/Avg_Q[n_tau][23][k])
		 <<"   "<<(Avg_Q[n_tau][18][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][19][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][20][k]/Avg_Q[n_tau][23][k])
		 <<"   "<<(Avg_Q[n_tau][21][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][22][k]/Avg_Q[n_tau][23][k])<<"   "<<Avg_Q[n_tau][23][k]<<endl;

		k=1;
		 out5<<fixed<<setprecision(10)<<Avg_Q[n_tau][0][k]<<"   "<<(Avg_Q[n_tau][1][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][2][k]/Avg_Q[n_tau][23][k])
         	<<"   "<<(Avg_Q[n_tau][3][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][4][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][5][k]/Avg_Q[n_tau][23][k])
		 <<"   "<<(Avg_Q[n_tau][6][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][7][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][8][k]/Avg_Q[n_tau][23][k])
		 <<"   "<<(Avg_Q[n_tau][9][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][10][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][11][k]/Avg_Q[n_tau][23][k])
		 <<"   "<<(Avg_Q[n_tau][12][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][13][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][14][k]/Avg_Q[n_tau][23][k])
		 <<"   "<<(Avg_Q[n_tau][15][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][16][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][17][k]/Avg_Q[n_tau][23][k])
		 <<"   "<<(Avg_Q[n_tau][18][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][19][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][20][k]/Avg_Q[n_tau][23][k])
		 <<"   "<<(Avg_Q[n_tau][21][k]/Avg_Q[n_tau][23][k])<<"   "<<(Avg_Q[n_tau][22][k]/Avg_Q[n_tau][23][k])<<"   "<<Avg_Q[n_tau][23][k]<<endl;	
	}	
	
}
