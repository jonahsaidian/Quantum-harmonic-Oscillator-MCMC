/**************************************************/
//    Atoms in a Magnetic Trap Simulation         //
//    For Physics 142 at UC San Diego             //
//    Code for Final Project                      //
/**************************************************/
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<numeric>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using std::cout;
using std::endl;

int main()
{
  const int N = 100;         //  Particle number
  int pos_bin[1000]= {0};     //  for binning location hits
  double dx = 0.01;           //  for good acceptance rate
  double k = 0.5;
  double T = .00010;      //  Includes Temp and Boltz Const
  double tempX[3];            //  Holds Test step
  srand((unsigned)time(NULL));
  FILE *rhist;
  rhist = fopen("./rhist.c","w");

  double **X;
  X = new double*[N];
  for(int i =0; i < N; i++)
    X[i] = new double[3];
  for(int i =0; i < N; i++)   //  Initialize Gas
    for(int j =0; j < 3; j++)
      X[i][j] =0.04*dx*((double)rand() / RAND_MAX - 0.5);

  long count = 0;
  long n = 0;
  while( n < 100000)
    {
      for (int i = 0; i < N; i++ )
	{
	  for(int j =0; j < 3; j++)   // step in 3D
	    tempX[j] = X[i][j] + dx*(2.0*(double)rand() / RAND_MAX -1);
	  double DE = k*(std::inner_product(tempX,tempX+3,tempX,0.0)
			 -std::inner_product(X[i],X[i]+3,X[i],0.0));
	  //  cout << DE << endl;   // for troubleshooting
	  //  cout << exp(-DE/T) << endl;
	  if ((double)rand()/RAND_MAX < exp(-DE/T))
	    for(int j =0; j < 3; j++)
	      X[i][j]=tempX[j];
	  if(X[i][1]== tempX[1])
	    count++;
	}
      n++;
      //cout << static_cast<double>(count)/(n*N)<< endl;  // for troubleshooting
    }
  count = 0;                        // re-clock
  long measure_count = 0;           // for average calculations
  double E,E2,Eave,Eerror,E2ave;    // for E and E2 ave on path and E ave over MC Sim
  double Etot = 0;                  // to hold E after every 100 steps
  double E2tot = 0;                 // to hold E sqrd for every 100 steps
  n =0;                             // reset n for convenience

  std::ofstream outFile;
  outFile.open("pos2.dat");
  outFile << std::setprecision(6) << std::fixed << std::showpoint;   // outputing waveform

  while( n < 1000000)
    {
      for (int i = 0; i < N; i++ )
	{
	  for(int j =0; j < 3; j++){
	    tempX[j] = X[i][j] + dx*(2.0*(double)rand() / RAND_MAX - 1.0);}
	  double DE = k*(std::inner_product(tempX,tempX+3,tempX,0.0)
			 -std::inner_product(X[i],X[i]+3,X[i],0.0));
	    //(std::inner_product(tempX, tempX+3, tempX, 0.0)
	    //		 -std::inner_product(X[i], X[i]+3, X[i], 0.0));
	  if ((double)rand()/RAND_MAX < exp(-DE/T))
	    for(int j =0; j < 3; j++)
	      X[i][j]=tempX[j];
	  if(X[i][1]== tempX[1])
	    count++;
	}
      if(0 == n %100)     // measure every 100
	{
	  measure_count++;
	  double Etemp;
	  E = 0;
	  E2 = 0;

	  for(int i = 0; i < N; i++)
	    {
	      double R = std::inner_product(X[i],X[i]+3,X[i],0.0);
	      Etemp = k*R;
	      R = sqrt(R);
	      E +=Etemp;
	      fprintf(rhist,"%f  %f  %f  \n",X[i][0],X[i][1],X[i][2]);
	      //      E2 += Etemp*Etemp;                   // average over gas
	      for(int j =0; j < 1000; j++)
		if(R >=j*.000003 && R < (j+1)*.000003)
		  pos_bin[j]++;
	    }
	  Etot += E;
	  E2tot += E*E;
	  Eave = Etot/(static_cast<double>(measure_count));               // Average over MC Sim
	  E2ave = E2tot/(static_cast<double>(measure_count));
	  Eerror = sqrt((E2ave - Eave*Eave)/(static_cast<double>(measure_count)));
	  if(0 == n%10000){
	    cout << "Percent done: " << n/10000.0 <<  "%" << endl;
	    cout << "Average energy: "<< Eave/N << "\tError in E is:  " << Eerror << endl;
	    cout << "Acceptance Rate:  " <<static_cast<double>(count)/(N*n) << endl;
	  }}
      n++;

    }
  for(int j =0; j < N; j++)
    delete X[j];
  delete X;
  cout << "\n\n\n" << endl;
  cout << "Analytical Solutions: " << 3.0*T/2.0 << endl;
  cout << "Finalaverage E is : " << Etot/(static_cast<double>(measure_count)*N) << endl;
  cout << "Final Error over Monte Carlo Simulation: " << Eerror << endl;
  for(int i = 0; i < 1000; i++)
    outFile <<  i*.000003 + 0.0000015 << "\t" << pos_bin[i] << endl;
  outFile.close();
  return 0;
  }
