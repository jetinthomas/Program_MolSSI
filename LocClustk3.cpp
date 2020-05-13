#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iomanip>

using namespace::std;

// Definining a maximum length to the arrays or matrices used in the code
#define T 1000
//Initializing Variables for the code, their role will come going further down the code
long int ind,N,M,ITE,Ite,L,Lold=10,point,sim;
static long double A[T],LnOme[T],H[T],Old[T],Aver[T][4],r,OrdParaDist[T][4*T],SuscDist[T][4*T],UnBlocked[T][5*T],Similar[T][4*T];
static int Lattice[T][T],B[T],Reference[T][4*T][T][T];
double Per;

// The void remember function stores the north-east/south-east diagonal before the particles in it is deleted. This is done to come back to the previous lattice if the new lattice formed is rejected in the metropolis algorithm.

void remember(int x, int y)
{
	// ind = 0 will store the north-east diagonal else south-east diagonal starting from 'x y' lattice site.

	if(ind==0){for(int i=0;i<=L-1;i++){B[i]=Lattice[x][y]; /* cout<<x<<"  "<<y<<"  "<<B[i]<<"  "<<endl;*/x=(L+x-1)%L; y=(y+1)%L;} /*cout<<endl;*/}

	else{for(int i=0;i<=L-1;i++){B[i]=Lattice[x][y]; /*cout<<x<<"  "<<y<<"  "<<B[i]<<"  "<<endl;*/x=(x+1)%L; y=(y+1)%L; } /*cout<<endl;*/}


}

// The void filling1 function fills the lattice site 'c i' and blocks the sites up to three nearest neighbors and increases the Number of Particles by 1.

void filling1(int c, int i)
{
				for(int d=-2;d<=2;d++)
							{for(int k=-2;k<=2;k++)
							{if(fabs(d)+fabs(k)<=2){Lattice[(L+c+d)%L][(L+i+k)%L]=Lattice[(L+c+d)%L][(L+i+k)%L]-1;}}
							}

				Lattice[c][i]=1; M++;
}

// The void deleting1 function deletes the lattice site 'c i' and unblocks the sites up to three nearest neighbors and decreases the Number of Particles by 1.

void deleting1(int c, int i)
{

							M=M-1;
							for(int d=-2;d<=2;d++)
							{for(int k=-2;k<=2;k++)
							{if(fabs(d)+fabs(k)<=2){Lattice[(L+c+d)%L][(L+i+k)%L]=Lattice[(L+c+d)%L][(L+i+k)%L]+1;}}
							}
							Lattice[c][i]=0;

}

// The void deleting function deletes all the particles in a north east diagonal if ind=0 or south-east diagonal otherwise.

void deleting(int i, int j)
{

	int l=L;
					if(ind==0){


                    //The loop runs till the Length of the north-east diagonal isn't traversed.
					while(l>0)
					{

						if(Lattice[i][j]==1)
						{
							M=M-1;
							for(int d=-2;d<=2;d++)
							{for(int k=-2;k<=2;k++)
							{if(fabs(d)+fabs(k)<=2){Lattice[(L+i+d)%L][(L+j+k)%L]=Lattice[(L+i+d)%L][(L+j+k)%L]+1;}}
							}
							Lattice[i][j]=0;i=(L+i-2)%L;j=(j+2)%L; l=l-2;
						}
						else{i=(L+i-1)%L;j=(j+1)%L;l--;}
					}
					}

					else{

                    //The loop runs till the Length of the south-east diagonal isn't traversed.
					while(l>0)
					{

						if(Lattice[i][j]==1)
						{
							M=M-1;
							for(int d=-2;d<=2;d++)
							{for(int k=-2;k<=2;k++)
							{if(fabs(d)+fabs(k)<=2){Lattice[(L+i+d)%L][(L+j+k)%L]=Lattice[(L+i+d)%L][(L+j+k)%L]+1;}}
							}
							Lattice[i][j]=0;i=(i+2)%L;j=(j+2)%L; l=l-2;
						}
						else{i=(i+1)%L;j=(j+1)%L;l--;}
					}}

}


// The void filling function fills particles in a unblocked segment in north east direction if ind=0 or south-east direction otherwise and blocks the sites up to third nearest neighbors from the filled sites.

void filling(int l, int c, int i)       //Inputs the length of the unblocked segment (l), and the site ('c i') where the segment begins.
{
	// For segment in 'North-East' direction
	if(ind==0){

	while(l>0)
	{

		if(l==L)
		{

			if(drand48()<= Per)     // 'Per' stores the probability for the site to get filled if the length of the segment is the lattice size, as due to periodic boundary condition this case considered specially.
			{

				for(int d=-2;d<=2;d++)
							{for(int k=-2;k<=2;k++)
							{if(fabs(d)+fabs(k)<=2){Lattice[(L+c+d)%L][(L+i+k)%L]=Lattice[(L+c+d)%L][(L+i+k)%L]-1;}}
							}

							Lattice[c][i]=1;c=(L+c-2)%L;i=(i+2)%L; l=l-3;
				 M++;
			}

			else
			{l=l-1; c=(L+c-1)%L; i=(i+1)%L;}
		}

		else
		{
			if(drand48()<= A[l])    //'A[l]' stores the probability for the site to get filled if length of the segment is 'l'.
			{
				for(int d=-2;d<=2;d++)
							{for(int k=-2;k<=2;k++)
							{if(fabs(d)+fabs(k)<=2){Lattice[(L+c+d)%L][(L+i+k)%L]=Lattice[(L+c+d)%L][(L+i+k)%L]-1;}}
							}
							Lattice[c][i]=1;c=(L+c-2)%L;i=(i+2)%L; l=l-2;
				M++;
            }
			else
			{l=l-1; i=(i+1)%L; c=(L+c-1)%L;}
		}

	}}

    // For segment in 'South-East' direction

	else{
		while(l>0)
	{
		if(l==L)
		{
			if(drand48() <= Per) // 'Per' stores the probability for the site to get filled if the length of the segment is the lattice size, as due to periodic boundary condition this case considered specially.
			{

				for(int d=-2;d<=2;d++)
							{for(int k=-2;k<=2;k++)
							{if(fabs(d)+fabs(k)<=2){Lattice[(L+c+d)%L][(L+i+k)%L]=Lattice[(L+c+d)%L][(L+i+k)%L]-1;}}
							}
							Lattice[c][i]=1;c=(c+2)%L;i=(i+2)%L; l=l-3;
				 M++;
			}

			else
			{l=l-1; c=(c+1)%L; i=(i+1)%L;}
		}

		else
		{
			if(drand48() <= A[l])//'A[l]' stores the probability for the site to get filled if length of the segment is 'l'.
			{

				for(int d=-2;d<=2;d++)
							{for(int k=-2;k<=2;k++)
							{if(fabs(d)+fabs(k)<=2){Lattice[(L+c+d)%L][(L+i+k)%L]=Lattice[(L+c+d)%L][(L+i+k)%L]-1;}}
							}
							Lattice[c][i]=1;c=(c+2)%L;i=(i+2)%L; l=l-2;
				M++;
            }
			else
			{l=l-1; i=(i+1)%L; c=(c+1)%L;}
		}

	}}

}

// The void comparing function rejects the new lattice whenever 'exp(LnOme[N]-LnOme[M])<z'

void comparing(int c, int i)
{

    double z=drand48();
	if((exp(1.0*(LnOme[N]-LnOme[M])))<z) // LnOme[] is an array stores density of states corresponding to each number of particle density. N is number of particles for old lattice and M is for new lattice.
	{

        // Modifying for north-east diagonal for ind=0 otherwise south-east diagonal. This modification is to get the old lattice back after the new lattice is rejected.
		if(ind==0)
		{
			deleting(c,i);  // Deletes the particles, the north-east diagonal starting from 'c i'. To, fill it again as per the old lattice.

		for(int i=0;i<=L-1;i++)
		{



			if(B[i]==1) // B[ ] is an array that stores the north-east diagonal starting from 'c i' of the old lattice in the void remember function.
			{
				for(int d=-2;d<=2;d++)
							{for(int k=-2;k<=2;k++)
							{if(fabs(d)+fabs(k)<=2){Lattice[(L+c+d)%L][(L+i+k)%L]=Lattice[(L+c+d)%L][(L+i+k)%L]-1;}}
							}
				Lattice[c][i]=1; M++;
			}
			c=(L+c-1)%L;
		}

		}

		else
		{

        // The same done for south-east diagonal if ind=1
		deleting(c,i);
		for(int c=0;c<=L-1;c++)
		{

			if(B[c]==1)
			{
				for(int d=-2;d<=2;d++)
							{
								for(int k=-2;k<=2;k++)
								{
									if(fabs(d)+fabs(k)<=2){Lattice[(L+c+d)%L][(L+i+k)%L]=Lattice[(L+c+d)%L][(L+i+k)%L]-1;}}
								}

				Lattice[c][i]=1; M++;
			}
			i=(i+1)%L;

		}

		}

	}

	if(point==1)   //The quantities are calculated when point=1 i.e. after the loop when histogram is within the required bound.
	{

		double NA[5],NB[5];
	    for(int fil=0;fil<5;fil++){NA[fil]=0;NB[fil]=0;}
	    for(int row=0;row<L;row++)
		{
			int sig=0;
			int col=(2*row)%5;
			int col1=(3+((3*row)%5))%5;
			while(sig<L)
			{

					if(Lattice[row][col]==1){NA[(5+col-((2*row)%5))%5]++;}

					if(Lattice[row][col1]==1){NB[(5+col1-(3+((3*row)%5))%5)%5]++;}

					sig++;


					col=(col+1)%L;	col1=(col1+1)%L;
			}
		}

		double qA=0,qB=0;
		for(int fil=0;fil<5;fil++)
		{
			int inc=fil+1;
			while(inc<5)
			{qA=qA+fabs(NA[fil]-NA[inc]);
			qB=qB+fabs(NB[fil]-NB[inc]);
			inc++;}
		}

		double q=(fabs(qA-qB));
		Aver[M][0]=Aver[M][0]+(1.0*M);               // For Particle Density
		Aver[M][1]=Aver[M][1]+pow((1.0*M),2);         // For Particle Density square
		Aver[M][2]=Aver[M][2]+q;              // For Order Parameter
		Aver[M][3]=Aver[M][3]+pow(q,2);       // For Order Parameter square


    if(OrdParaDist[M][int(q)]==0)       //Storing the reference lattice for a particular density and packing fraction
    {

      // cout<<"Reference Lattice for M and q "<<OrdParaDist[M][int(q)]<<"  "<<M<<"   "<<int(q)<<"  "<<Similar[M][int(q)]<<endl;

        for(int i=0;i<L;i++)

            {
                for(int j=0; j<L; j++)
                   {
                       Reference[M][int(q)][i][j]=Lattice[i][j];
                       // cout<<Reference[M][int(q)][i][j]<<"  ";
                   }
               // cout<<endl;
            }

       // cout<<endl<<endl;
    }

    else
        {
            //Similar[M][int(q)]=0;
            for(int i=0;i<L;i++)                //Calculating how many same sites are occupied  for reference and constructed lattice.

            {
                for(int j=0; j<L; j++)
                   {
                       if(Reference[M][int(q)][i][j]==1 && Lattice[i][j]==1)
                       {Similar[M][int(q)] = Similar[M][int(q)]+ 1;}

                   }
            }


        }


	OrdParaDist[M][int(q)] = 1.0 + OrdParaDist[M][int(q)];  // For OrderParameter Distribution

	int cont=0;
	for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++){if(Lattice[i][j]==0){cont++;}}
    }

    UnBlocked[M][cont] = UnBlocked[M][cont]+1;    //For UnBlocked Sites Distribution

	}
		H[M]=H[M]+1;    //Updating the histogram once we get the required lattice after a acceptance-rejection move.

}


void calculate()
{

    double NA[5],NB[5];
	for(int fil=0;fil<5;fil++){NA[fil]=0;NB[fil]=0;}
	for(int row=0;row<L;row++)
	{
		int sig=0;
		int col=(2*row)%5;
		int col1=(3+((3*row)%5))%5;
			while(sig<L)
			{

					if(Lattice[row][col]==1){NA[(5+col-((2*row)%5))%5]++;}

					if(Lattice[row][col1]==1){NB[(5+col1-(3+((3*row)%5))%5)%5]++;}

					sig++;

										col=(col+1)%L;	col1=(col1+1)%L;
			}
	}

	double qA=0,qB=0;
	for(int fil=0;fil<5;fil++)
	{
		int inc=fil+1;
		while(inc<5)
		{qA=qA+fabs(NA[fil]-NA[inc]);
		qB=qB+fabs(NB[fil]-NB[inc]);
		inc++;}
	}

	double q=fabs(qA-qB);

	Aver[M][0]=Aver[M][0]+(1.0*M);        // For Particle Density
	Aver[M][1]=Aver[M][1]+pow((1.0*M),2); // For Particle Density square
	Aver[M][2]=Aver[M][2]+q;      // For Order Parameter
	Aver[M][3]=Aver[M][3]+pow(q,2);   // For Order Parameter square


    if(OrdParaDist[M][int(q)]==0)       //Storing the reference lattice for a particular density and packing fraction
    {
         //cout<<"Reference Lattice for M and q "<<OrdParaDist[M][int(q)]<<"  "<<M<<"   "<<int(q)<<endl;

            for(int i=0;i<L;i++)

            {
                for(int j=0; j<L; j++)
                   {
                       Reference[M][int(q)][i][j]=Lattice[i][j];
                       // cout<<Reference[M][int(q)][i][j]<<"  ";
                   }
               // cout<<endl;
            }
      //  cout<<endl;
    }

    else
        {

           // Similar[M][int(q)]=0;
            for(int i=0;i<L;i++)                    //Calculating how many same sites are occupied  for reference and constructed lattice.
            {
                for(int j=0; j<L; j++)
                   {
                       if(Reference[M][int(q)][i][j]==1 && Lattice[i][j]==1)
                       {Similar[M][int(q)] = Similar[M][int(q)]+ 1;}

                   }
            }


        }

	OrdParaDist[M][int(q)] = 1.0+OrdParaDist[M][int(q)];  // For Order Parameter Distribution


	int cont=0;
	for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++){if(Lattice[i][j]==0){cont++;}}
    }

    UnBlocked[M][cont] = UnBlocked[M][cont]+1;    // For UnBlocked Sites Distribution
}



int main()
{
	srand48(152211);
	string Result;
	string Result1;

	long int Size[6];
	Size[0]=10;Size[1]=15;Size[2]=20;Size[3]=25;Size[4]=30;Size[5]=40;      // Definining the Length of lattices

	for(int jump=0;jump<=5;jump++)      // loop for different sizes of lattices
	{
	    L=Size[jump];
		ostringstream convert;
		convert<<L;
		Result = convert.str();
		ifstream inp ("obc_probk3L1000mu0.0000",ios::out);        // Reading file that stored probabilities for open segments
		ifstream inp1 ("pbc_probk3L1000mu0.0000",ios::out);       // Reading file that stored probabilities for periodic segments

		for(int i=0;i<=L;i++)
		{
			inp>>i>>A[i];               //Storing prob for open segments into A[]
			inp1>>i>> Per ;             //Storing prob for periodic segments into Per
		}

		long double Hbar;
    	double prob;
    	if(L<=20){prob=0.1;}    // Defining probabilities for single site moves.
    	else{prob=0.4;}

    	ostringstream convert1;
    	convert1<<prob;
    	Result1=convert1.str();

		ofstream out(("L"+Result+"p"+Result1+"_Evolution_2x10^(6).dat").c_str(),ios::out);      // Opens a file for storing histogram after each iteration
		ofstream out1(("L"+Result+"p"+Result1+"_Averages_2x10^(6).dat").c_str(),ios::out);      // Open a file for storing average quantities after the final iteration
    	ofstream out2(("L"+Result+"p"+Result1+"_OrdParaDist_2x10^(6).dat").c_str(),ios::out);   // Open a file for storing Order Parameter Distribution after the final iteration
    	ofstream out3(("L"+Result+"p"+Result1+"_UnBlockSitesDist_2x10^(6).dat").c_str(),ios::out); // Open a file for storing UnBlocked Site Distribution after the final iteration
    	ofstream out4(("L"+Result+"p"+Result1+"_SimilarConfDist_2x10^(6).dat").c_str(),ios::out); // Open a file for storing UnBlocked Site Distribution after the final iteration

  	   // Headings for all the opened files

    	out<<"Density"<<"  "<<"Normalised"<<"  "<<"Iteration"<<"  "<<"Fraction of"<<"  "<<"Ln(Density"<<endl;
    	out<<"       "<<"  "<<"Histogram"<<"  "<<"         "<<"  "<<"Accepted Moves"<<"  "<<"of states)"<<endl<<"#"<<endl;
    	out1<<"Density"<<"  "<<"Normalised"<<"  "<<"Ln(Density"<<"  "<<"Average"<<"  "<<"Average"<<"  "<<"Average"<<"   "<<"Average"<<endl;
    	out1<<"       "<<"  "<<"Histogram"<<"  "<<" of states)"<<"  "<<"Density"<<"  "<<"Density^2"<<"  "<<"q"<<"   "<<"q^2"<<endl<<"#"<<endl;
    	out2<<"Density"<<"  "<<"OrderParameter"<<"  "<<"OrderParameter Dist"<<endl<<"#"<<endl;
    	out3<<"Density"<<"  "<<"UnBlock Sites"<<"  "<<"UnBlock Sites Dist"<<endl<<"#"<<endl;
    	out4<<"Density"<<"  "<<"Order Parameter"<<"  "<<"Similar Conf Dist"<<endl<<"#"<<endl;

    	// Initialising Matrices that will store the density, order-parameter and the distributions.

    	for(int i=0;i<=(L*L/5);i++)
		{
	    	Aver[i][0]=0;Aver[i][1]=0;Aver[i][2]=0;Aver[i][3]=0;
	    	for(int j=0;j<=(L*L);j++)
			{
        		if(j<=4*L*L/5){OrdParaDist[i][j]=0; SuscDist[i][j]=0;Similar[i][j]=0;} UnBlocked[i][j]=0;
			}
		}

		point=0;
		Ite=1;  // The Iteration
    int test=0;
	while(point<2)  // The loop that will run till the histogram is within the required bound.
	{
        ITE=0;
		double POW=2*pow(10,6); // The base number of steps in every iteration
		int SIM=10;              // The number of initial configurations
		for(int i=0;i<=0.20*L*L;i++)
		{
		    // The loop that updates Density of states 'LnOme[]'  and histogram 'H[]' after each iteration ('Ite') but not for the first iteration.
			if(Ite>1){LnOme[i]=LnOme[i]+log(1.0*(H[i]+1)/(Hbar/((0.20*L*L+1)))); H[i]=0;}

            // For Ite =1
			else{

                    // For other size than the minimum, which initializes the density of states from the one obtained from smaller lattice.
					if(jump!=0)
					{int den2=0;
					for(int den1=0;den1<=(L*L)/5;den1++)
					{
						while(1.0*den2/(Lold*Lold)<1.0*den1/(L*L)){den2++;}

						if(den1!=0){LnOme[den1] = ((1.0*L*L)/(Lold*Lold))*(Old[den2-1]+(Lold*Lold)*(Old[den2]-Old[den2-1])*((1.0*den1/(L*L))-((1.0*den2-1)/(Lold*Lold))));}
						else{LnOme[den1] = ((1.0*L*L)/(Lold*Lold))*Old[den2];}

					}}

                    // For the minimum size lattice.
					else{LnOme[i]=0;}

					H[i]=0; // Histogram is initialized to zero.
				}

		}

		Hbar=0; // Hbar counts the total number of samples in a iteration, this is also initialized to zero.

		for(sim=0;sim<=SIM-1;sim++){

        // Initializing the Lattice before the iteration begins.
		for(int i=0;i<L;i++)
		{
			for(int j=0;j<L;j++)
			{
				Lattice[i][j]=0;

			}
		}


		M=0; // Initializing for the 'M' for the particle number




    // The initial lattice is created depending upon 'sim'. The particle density will vary from 0 to 0.2 based on sim.
	for(int row=0;row<L;row++)
	{
        int col;
		if(M<(sim%5)*L*L/20)
        {
                if(sim<5){col=(2*row)%5;}
                else{col=(3+((3*row)%5))%5;}
				for(int d=-2;d<=2;d++)
                {for(int k=-2;k<=2;k++)
                {if(fabs(d)+fabs(k)<=2){Lattice[(L+row+d)%L][(L+col+k)%L]=Lattice[(L+row+d)%L][(L+col+k)%L]-1;}}}
                Lattice[row][col]=1;
                M++;
        }
	}

		N=0;

		r=0;

        // The iteration begins

		while(r<2*Ite*POW)
		{

            // prob decides to do the single site move or updating a diagonal.
			if(drand48()<prob){

					int c,i;
					// Selects a site randomly.
					c=int(1.0*(L)*drand48());
					i=int(1.0*(L)*drand48());
					int skip=0;
                            // Fills, Deletes or Keep as it is if the chosen site is blocked.
							if(Lattice[c][i]==1)
									{

										if(drand48()<exp(LnOme[M]-LnOme[M-1])){deleting1(c,i);}  skip++; // deletes

									}

									if(Lattice[c][i]==0 && skip==0)
									{
										if(drand48()<exp(LnOme[M]-LnOme[M+1])){filling1(c,i);} // fills

									}

										Hbar=Hbar+1;H[M]++;r++; // Updating
										if(point==1){calculate();} // Going to function calculate() to get the average quantities
                                }
				else{
					int c,i;
					// Selects a diagonal randomly, either it would be traversing north-east or south-east
					if(drand48()<=0.5){c=int(1.0*(L)*drand48());i=0;ind=0;}
					else{c=0;i=int(1.0*(L)*drand48());ind=1;}

                    // For North-East
					if(ind==0)
					{
								N=M;
								remember(c,i); // Stores the diagonal
								deleting(c,i);  // Deletes the particles in the diagonal
								int l=0,skip=0;
								int p=c;
								while(i<L)
								{
								    // If the entire diagonal is unblocked
									while(Lattice[c][i]==0 && i<L)
									{l++;c=(L+c-1)%L;i=(i+1);skip=1;}

                                     // Identifies segments in the diagonal having no blocked sites and their site in-dices
									if(l!=L)
									{
										l=0;
										while(Lattice[c][i]!=0 && i<L)
										{c=(L+c-1)%L;i=(i+1);}
										int ti1=c%L,ti2=i%L;
										while(Lattice[ti1][ti2]==0)
										{ti1=(L+ti1-1)%L;ti2=(ti2+1)%L;l++;skip=1;}

									}

									filling(l,c,i%L); // fills the segment
									i=i+l;c=(L+c-l)%L;

                                }

							comparing(p,0); // Calls the comparing function for accepting or rejecting the new lattice formed.
							Hbar++; // Updating

						}

                        // For South-east

						else
						{

								N=M;
								remember(c,i); // Stores the diagonal
								deleting(c,i); // Deletes the particles in the diagonal
								int l=0,skip=0;
								int p=i;
								while(c<L)
								{
								    // If the entire diagonal is unblocked
									while(Lattice[c][i]==0 && c<L)
									{l++;c=(c+1);i=(i+1)%L;skip=1;}

                                    // Identifies segments in the diagonal having no blocked sites and their site in-dices
									if(l!=L)
									{
										l=0;
										while(Lattice[c][i]!=0 && c<L)
										{c=c+1;i=(i+1)%L;}
										 int ti1=c%L,ti2=i%L;
										while(Lattice[ti1][ti2]==0)
										{ti1=(ti1+1)%L;ti2=(ti2+1)%L;l++;skip=1;}

									}

									filling(l,c%L,i); // fills the segment
									i=(i+l)%L;c=c+l;

							}
							comparing(0,p); // Calls the comparing function for accepting or rejecting the new lattice formed.
							Hbar++;  // Updating
						}

					r=r+1;

					}

		}
}

    // Output the histogram at different iterations
	for(int i=0;i<=0.20*L*L;i++){out<<1.0*i/(L*L)<<"  "<<(H[i])/(Hbar/((0.20*L*L+1)))<<"  "<<Ite<<"  "<<(1.0*Hbar/(r*SIM))<<"  "<<(LnOme[i]-LnOme[0])<<endl;}out<<"#"<<endl<<"#"<<endl;
    // Estimates the histogram is within a bound.
    for(int i=0;i<=0.20*L*L;i++){if(0.8<(H[i])/(Hbar/((0.20*L*L+1))) && (H[i])/(Hbar/((0.20*L*L+1))) < 1.2){ITE++;}}
    //Estimate to stop the loop after one more run after this or to continue for more loops
    if(ITE==(L*L)/5+1 || point==1){point++;}
	else{Ite++;}
	if(point==2){
        // Output for final loop
        for(int i=0;i<=0.20*L*L;i++){out1<<1.0*i/(L*L)<<"  "<<(H[i])/(Hbar/((0.20*L*L+1)))<<"  "<<(LnOme[i]-LnOme[0])<<"   "<<pow((1.0/(L*L)),1)*(Aver[i][0]/H[i])
        <<"  "<<pow((1.0/(L*L)),2)*(Aver[i][1]/H[i])<<"   "<<pow((5.0/(4*L*L)),1)*(Aver[i][2]/H[i])<<"  "<<pow((5.0/(4*L*L)),2)*(Aver[i][3]/H[i])<<endl;}



    for(int i=0;i<=L*L/5;i++){double sum=0;
    for(int j=0;j<=4*L*L/5;j++)
    {sum=sum+OrdParaDist[i][j];if(OrdParaDist[i][j]!=0){out2<<1.0*i/(L*L)<<"    "<<5.0*j/(4*L*L)<<"    "<<OrdParaDist[i][j]<<endl;}}
        out2<<"#"<<endl<<"#"<<endl;out2<<"Normalised Cumalitive Sum"<<"   "<<sum<<"   "<<H[i]<<endl<<"#"<<endl;}

    for(int i=0;i<=L*L/5;i++){double sum=0;
    for(int j=0;j<=L*L;j++)
        {sum=sum+UnBlocked[i][j];if(UnBlocked[i][j]!=0){out3<<1.0*i/(L*L)<<"    "<<j<<"    "<<UnBlocked[i][j]<<endl;}}
        out3<<"#"<<endl<<"#"<<endl;out3<<"Normalised Cumalitive Sum"<<"   "<<sum<<"   "<<H[i]<<endl<<"#"<<endl;}

    for(int i=0;i<=L*L/5;i++){
    for(int j=0;j<=4*L*L/5;j++)
    {
		if(OrdParaDist[i][j]!=0 && OrdParaDist[i][j]!=1){out4<<1.0*i/(L*L)<<"    "<<5.0*j/(4*L*L)<<"    "<<1.0*Similar[i][j]/(L*L*(OrdParaDist[i][j]-1))<<endl;}}
        }

        }
	}
    // Stored the estimated density of states for using it for giving the initial estimate of density of states for the next larger lattice.
	for(int i=0;i<=0.20*L*L;i++){Old[i]=LnOme[i]-LnOme[0];}
	Lold=L;
}
}
