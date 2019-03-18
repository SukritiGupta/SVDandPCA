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
		cout<<*(m+i*N+j);
	}
	// cout<<endl<<len<<endl;

	len=sqrt(len);
	*(r+i*N+i)=len;
	cout<<endl<<"make unit row"<<len<<endl;
	if(len!=0)
		for (int j = 0; j < N; ++j)
		{
			*(m+i*N+j)=*(m+i*N+j)/len;
		}

}
double getdotc(double*a, int j, double*q, int i)
{
	double ret=0.0;
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
void QRFactors(double*a, double*q, double*r) ///write this again!!!!!!!!!!!!!1
{
	using namespace std;
	for (int i = 0; i < N; ++i)
	{
		makeunitrow(q,i,r);
		for (int j = i+1; j < N; ++j)
		{
			*(r+i*N+j)=getdot(a,j,q,i);
			
		}

		for (int i1 = i+1; i1 < N; ++i1)
		{
		// int i1=1;
			double valdot=*(r+i*N+i1);
			for (int j1 = 0; j1 < N; ++j1)
			{
				// cout<<endl<<(*(q+i1*N+j1))<<endl<<valdot<<endl<<*(q+i*N+j1)<<endl;
				(*(q+i1*N+j1))=(*(q+i1*N+j1))-(valdot)*(*(q+i*N+j1));
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
	using namespace std;
	double* dtd;
	dtd = (double*) malloc(sizeof(double) * N*N);	
	int p=1;

	for (int i = 0; i < N; i+=1)
	{
		int temp=i*N;
		for (int j = i; j < N; ++j)
		{
			*(dtd+temp+j)=0;			
		}		
	}

	// 1. strided access-d 2. Cache line invalidation - dtd
	for (int o = 0; o < M; ++o)
	{
		for (int i = 0; i < N; i+=1)
		{
			int temp=i*N;
			for (int j = i; j < N; ++j)
			{
				*(dtd+temp+j)+=(*(D+o*N+i))*(*(D+o*N+j));	
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

	for (int i = 0; i < N; i+=1)
	{
		int temp=i*N;
		for (int j = 0; j < N; ++j)
		{
			*(r+temp+j)=0;
			*(e+temp+j)=0;
		}		
	}

	for (int i = 0; i < N; ++i)
	{
		*(e+i*N+i)=1;
	}



	int loop=0,flag=0;
	while(true)
	{

		for (int i = 0; i < N; i+=1)
		{
			int temp=i*N;
			for (int j = 0; j < N; ++j)
			{
				*(q+i+j*N)=*(dtd+temp+j);			
			}		
		}

		QRFactors(dtd, q,r);

		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				*(dtd+i*N+j)=getdotc(r,i,q,j);
				*(e1+i*N+j)=getdotc(e,i,q,j);
				if(flag==0 && i!=j)
				{
					if (*(dtd+i*N+j)>0.00001||*(dtd+i*N+j)<-0.00001)
					{
						flag=1;
					}
				}	
			}
		}

		// e=e1;

		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				*(e+i*N+j)=*(e1+i*N+j);
			}
		}
		loop++;
		if(flag==0||loop>50)
		{
			break;
		}

	}


	cout<<endl<<"eigenval"<<endl;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			printf("%f\n", *(dtd+i*N+j));
		}
		printf("\n");
	}

	cout<<endl<<"eigenvector"<<endl;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			printf("%f\n", *(e+i*N+j));
		}
		printf("\n");
	}

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
  		cout<<(*(rit)).lam<<endl;
  		cout<<(*(rit)).i<<endl;
  		*(*SIGMA+ii)=sqrt((*(rit)).lam);

  		ii++;

  	}

  	for (int i = 0; i < N; ++i)
  	{
  		// *(*SIGMA+i)=g1.at(i).lam;
  		// cout<<*(*SIGMA+i);
  		printf("%f\n",*(*SIGMA+i) );
  	}
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
  	for (int i = 0; i < M; ++i)
  	{
  		for (int j = 0; j < M; ++j)
  		{
  			*(merau+i*M+j)=0;
  		}
  	}


  	for (int i = 0; i < M; ++i)
  	{
  		for (int j = 0; j < N; ++j)
  		{
  			if ((*(*SIGMA+j))!=0)
  				*(merau+i*M+j)=(getmul(D,i,*U,j))/(*(*SIGMA+j));
  			cout<<"assigned"<<i<<j<<endl;
  		}
  		cout<<"output merau"<<endl;
	  	for (int i1 = 0; i1 < M; ++i1)
	  	{
	  		for (int j1 = 0; j1 < M; ++j1)
	  		{
	  			cout<<*(merau+i1*M+j1)<<endl;
	  		}
	  	}





  	}

 //*******************************************************************


  	for (int i = 0; i < M; ++i)
  	{
  		for (int j = 0; j < M; ++j)
  		{
  			*(*V_T+i*M+j)=0;
  		}
  	}

  	  	cout<<"output V_T"<<endl;
  	for (int i = 0; i < M; ++i)
  	{
  		for (int j = 0; j < M; ++j)
  		{
  			cout<<*(*V_T+i*M+j)<<endl;
  		}
  	}

  	for (int i = 0; i < M; ++i)
  	{
  		for (int j = 0; j < M; ++j)
  		{
  			// float temp=(*(*SIGMA+j));
  			// if (temp!=0)
  			// {
  			// 	*(*V_T+j*N+i)=(getmul(D,i,*U,j))/temp;
  				
  			// }
  			*(*V_T+j*M+i)=*(merau+i*M+j);
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


  	cout<<"output U"<<endl;
  	for (int i = 0; i < N; ++i)
  	{
  		for (int j = 0; j < N; ++j)
  		{
  			cout<<*(*U+i*N+j)<<endl;
  		}
  	}

  	cout<<"output SIGMA"<<endl;
  	for (int i = 0; i < N; ++i)
  	{
  			cout<<*(*SIGMA+i)<<endl;
  	}

  	cout<<"output V_T"<<endl;
  	for (int i = 0; i < M; ++i)
  	{
  		for (int j = 0; j < M; ++j)
  		{
  			cout<<*(*V_T+i*M+j)<<endl;
  		}
  	}
}
