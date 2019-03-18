//correct 

#include <malloc.h>
#include <omp.h>
#include<iostream>
#include <math.h>
#include <vector>
#include <algorithm>

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
int N;
struct indlam
{
	double lam;
	int i;	
};

bool my_cmp(const indlam& a, const indlam& b)
{
    // smallest comes first
    return a.lam < b.lam;
}
float getmul(float*a, int i, float*b, int j)
{
	double ret=0.0;
	for (int k= 0; k < N; ++k)
	{
		ret+=*(a+i*N+k)*(*(b+k*N+j));
	}
	return (float)ret;
}
void makeunitrow(double*m, int i, double*r)
{
	using namespace std;
	double len=0;
	for (int j = 0; j < N; ++j)
	{
		len+=*(m+i*N+j)*(*(m+i*N+j));
		// cout<<*(m+i*N+j);
	}
	// cout<<endl<<len<<endl;

	len=sqrt(len);
	*(r+i*N+i)=len;
	// cout<<endl<<"make unit row"<<len<<endl;
	if(len!=0)
		for (int j = 0; j < N; ++j)
		{
			*(m+i*N+j)=*(m+i*N+j)/len;
		}

}
double getdotc(double*a, int j, double*q, int i)
{
	double ret=0.0;
	// #pragma omp parallel for reduction(+:ret)
	for (int k= 0; k < N; ++k)
	{
		ret+=*(a+j*N+k)*(*(q+i*N+k));
	}
	return (double)ret;
}
double getdot(double*a, int j, double*q, int i)
{
	double ret=0.0;
	for (int k= 0; k < N; ++k)
	{
		ret+=*(a+j+k*N)*(*(q+i*N+k));
	}
	return (double)ret;
}
void QRFactors(double*a, double*q, double*r) 
{
	using namespace std;
	for (int i = 0; i < N; ++i)
	{
		makeunitrow(q,i,r);
		for (int j = i+1; j < N; ++j)
		{
			*(r+i*N+j)=getdot(a,j,q,i);
			
		}
		#pragma omp parallel for collapse(2) 
		for (int i1 = i+1; i1 < N; ++i1)
		{
		// int i1=1;
			// double valdot=*(r+i*N+i1);
			for (int j1 = 0; j1 < N; ++j1)
			{
				// cout<<endl<<(*(q+i1*N+j1))<<endl<<valdot<<endl<<*(q+i*N+j1)<<endl;
				(*(q+i1*N+j1))=(*(q+i1*N+j1))-(*(r+i*N+i1))*(*(q+i*N+j1));
				// exit(0);
			}
			
		}
		//set value of e_i+1
	}
}




int Cache=16;
void SVD(int M, int gN, float* D, float** U, float** SIGMA, float** V_T)
{
	N=gN;
	int omax=M*M/20000;
	int max=5000;
	if(max<omax)
	{
		max=omax;
	}
	using namespace std;
	double* dtd;
	dtd = (double*) malloc(sizeof(double) * N*N);	
	int p=1;
	#pragma omp parallel for collapse(2)
	for (int i = 0; i < N; i+=1)
	{
		// int temp=i*N;
		for (int j = 0; j < N; ++j)
		{
			*(dtd+i*N+j)=0;			
		}		
	}

	// 1. strided access-d 2. Cache line invalidation - dtd
	// #pragma omp parallel for 
	for (int o = 0; o < M; ++o)
	{
		for (int i = 0; i < N; i+=1)
		{
			// int temp=i*N;
			for (int j = i; j < N; ++j)
			{
				*(dtd+i*N+j)+=(*(D+o*N+i))*(*(D+o*N+j));	
			}		
		}
	}

	for (int i = 1; i < N; ++i)
	{
		for (int j = i-1; j >=0; --j)
		{
			*(dtd+i*N+j)=*(dtd+N*j+i);
		}
	}


	double *q;
	double *r;
	double *e;
	double *e1;

	float *merau;


	q = (double*) malloc(sizeof(double) * N*N);	
	r = (double*) malloc(sizeof(double) * N*N);
	e = (double*) malloc(sizeof(double) * N*N);
	e1 = (double*) malloc(sizeof(double) * N*N);
	merau = (float*) malloc(sizeof(float) * M*M);

	#pragma omp parallel for collapse(2)
	for (int i = 0; i < N; i+=1)
	{
		// int temp=i*N;
		for (int j = 0; j < N; ++j)
		{
			*(r+i*N+j)=0;
			*(e+i*N+j)=0;
		}		
	}

	for (int i = 0; i < N; ++i)
	{
		*(e+i*N+i)=1;
	}



	int loop=0;
	
	while(true)
	{
		bool flag=false;

		#pragma omp parallel for collapse(2)
		for (int i = 0; i < N; i+=1)
		{
			// int temp=i*N;
			for (int j = 0; j < N; ++j)
			{
				*(q+i+j*N)=*(dtd+i*N+j);			
			}		
		}

		QRFactors(dtd, q,r);
		#pragma omp parallel
		{ 
			// cout<<omp_get_num_threads();

		#pragma omp for collapse(2) reduction(||:flag)
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				*(dtd+i*N+j)=getdotc(r,i,q,j);
				*(e1+i*N+j)=getdotc(e,i,q,j);
				if(flag==false && i!=j)
				{
					if (*(dtd+i*N+j)>0.0001||*(dtd+i*N+j)<-0.0001)
					{
						flag=true;
					}
				}	
			}
		}

		// e=e1;
		#pragma omp for collapse(2)
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				*(e+i*N+j)=*(e1+i*N+j);
			}
		}

		}
		loop++;
		
		// cout<<loop<<endl;
		if(flag==false||loop>max)
		{
			break;
		}

	}


	// cout<<endl<<"eigenval"<<endl;
	// for (int i = 0; i < N; ++i)
	// {
	// 	for (int j = 0; j < N; ++j)
	// 	{
	// 		printf("%f\n", *(dtd+i*N+j));
	// 	}
	// 	printf("\n");
	// }

	// cout<<endl<<"eigenvector"<<endl;
	// for (int i = 0; i < N; ++i)
	// {
	// 	for (int j = 0; j < N; ++j)
	// 	{
	// 		printf("%f\n", *(e+i*N+j));
	// 	}
	// 	printf("\n");
	// }

	vector<indlam> g1; 
  

	for (int i = 0; i < N; ++i)
	{
		indlam t={*(dtd+i*N+i),i};
		g1.push_back(t);		
	}
	
	sort(g1.begin(),g1.end(),my_cmp);
	//setting sigma correctly
	std::vector<indlam>::reverse_iterator rit = g1.rbegin();
	int ii=0;
  	for (; rit!= g1.rend(); ++rit)
  	{
  		// cout<<(*(rit)).lam<<endl;
  		// cout<<(*(rit)).i<<endl;
  		*(*SIGMA+ii)=sqrt((*(rit)).lam);

  		ii++;

  	}

  	// for (int i = 0; i < N; ++i)
  	// {
  	// 	// *(*SIGMA+i)=g1.at(i).lam;
  	// 	// cout<<*(*SIGMA+i);
  	// 	printf("%f\n",*(*SIGMA+i) );
  	// }
    // *rit = ++i;


  	//setting inka U correctly
    std::vector<indlam>::reverse_iterator it = g1.rbegin();
	int iii=0;
  	for (; it!= g1.rend(); ++it)
  	{
  		int index=(*it).i;
  		for (int i = 0; i < N; ++i)
  		{
  			*(*U+i*N+iii)=(float)(*(e+i*N+index));  			
  		}
  		iii++;
  	}

//***************************************************************
#pragma omp parallel
{
  	#pragma omp for collapse(2)
  	for (int i = 0; i < M; ++i)
  	{
  		for (int j = 0; j < M; ++j)
  		{
  			*(merau+i*M+j)=0;
  		}
  	}

  	#pragma omp for collapse(2)
  	for (int i = 0; i < M; ++i)
  	{
  		for (int j = 0; j < N; ++j)
  		{
  			if ((*(*SIGMA+j))!=0)
  				*(merau+i*M+j)=(getmul(D,i,*U,j))/(*(*SIGMA+j));
  			// cout<<"assigned"<<i<<j<<endl;
  		}
  	}

 //*******************************************************************

  	#pragma omp for collapse(2)
  	for (int i = 0; i < M; ++i)
  	{
  		for (int j = 0; j < M; ++j)
  		{
  			*(*V_T+i*M+j)=0;
  		}
  	}

  	#pragma omp for collapse(2)
  	for (int i = 0; i < M; ++i)
  	{
  		for (int j = 0; j < M; ++j)
  		{
  			*(*V_T+j*M+i)=*(merau+i*M+j);
  		}
  	}
}




  	//**************************************************************************************final check

  	// cout<<"input D"<<endl;
  	// for (int i = 0; i < M; ++i)
  	// {
  	// 	for (int j = 0; j < N; ++j)
  	// 	{
  	// 		cout<<*(D+i*N+j)<<endl;
  	// 	}
  	// }


  	// cout<<"output U"<<endl;
  	// for (int i = 0; i < N; ++i)
  	// {
  	// 	for (int j = 0; j < N; ++j)
  	// 	{
  	// 		cout<<*(*U+i*N+j)<<endl;
  	// 	}
  	// }

  	// cout<<"output SIGMA"<<endl;
  	// for (int i = 0; i < N; ++i)
  	// {
  	// 		cout<<*(*SIGMA+i)<<endl;
  	// }

  	// cout<<"output V_T"<<endl;
  	// for (int i = 0; i < M; ++i)
  	// {
  	// 	for (int j = 0; j < M; ++j)
  	// 	{
  	// 		cout<<*(*V_T+i*M+j)<<endl;
  	// 	}
  	// }
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
	using namespace std;
	float*W;
	int k=0;
	float sum=0;
	for (int i = 0; i < N; ++i)
	{
		sum+=(*(SIGMA+i))*(*(SIGMA+i));
		// cout<<*(SIGMA+i);
		// k++;
		// if (sum>=retention)
		// {
		// 	*K=k;
		// 	break;
		// }
	}

	float temp=0;
	for (int i = 0; i < N; ++i)
	{
		temp+=(*(SIGMA+i))*(*(SIGMA+i))/sum;
		// cout<<*(SIGMA+i);
		k++;
		if ((temp*100)>=retention)
		{
			*K=k;
			break;
		}
	}
	W = (float*) malloc(sizeof(float) * N*k);
	*(D_HAT)=(float*) malloc(sizeof(float) * M*k);

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < k; ++j)
		{
			*(W+i*k+j)=*(U+i*N+j);
		}
	}

	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < k; ++j)
		{
			float assign=0;
			for (int kk= 0; kk < N; ++kk)
			{
				assign+=(*(D+i*N+kk))*(*(W+kk*k+j));
			}
			*(*D_HAT+i*k+j)=assign;
		}
	}

	// cout<<k<<endl;



	// cout<<"W"<<endl;
	// for (int i = 0; i < N; ++i)
	// {
	// 	for (int j = 0; j < k; ++j)
	// 	{
	// 		cout<<*(W+i*k+j);
	// 	}
	// }

	// cout<<"D_HAT"<<endl;
	// for (int i = 0; i < 10; ++i)
	// {
	// 	for (int j = 0; j < 10; ++j)
	// 	{
	// 		cout<<*(*D_HAT+i*k+j)<<endl;
	// 	}
	// 	cout<<endl;
	// }
    
}

//complete sequential code on 1000*50 dataset took 20s
