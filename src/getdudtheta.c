#include "header.h" 


//START---getdmudtheta_ACD----------------------------//
void getdmudtheta_ACD(double *x,
		int *N,
		double *par,
		int *order,
		double *mean,
		double *mu,
		double *resi,
		int *newDay,
		int *NnewDays,
		double *dmudtheta)
{

	int maxpq = max(*order,*(order+1));
	int i = 0;
	int j = 1;
	int k = 0;
	int nextND=0;

	int startIndex = 0, stopIndex = maxpq;

	do {
		//in the start of the sample or in the start of a new day, the maxpq mus are set to the mean:
		for (i = startIndex; i < stopIndex; i++){
			*(mu+i) = *mean;
			*(resi+i) = *(x+i)/(*(mu+i));
			for(j = 0; j < 1+ *order + *(order+1); j++){
				*(dmudtheta + j * *N + i) = 0;
			}
		}
		startIndex = stopIndex;
		if(nextND < *NnewDays) stopIndex = *(newDay + nextND++) - 1;
		else stopIndex = *N;
		for(i = startIndex; i < stopIndex; i++){
			*(mu+i) = *par; //adds the constant
			for(j = 1; j <= *order; j++) *(mu+i) += *(par+j)*(*(x+i-j)); //adds the lagged durations
			for(j = 1; j <= *(order+1); j++) *(mu+i) += *(par+*order+j)*(*(mu+i-j)); //adds the lagged mus
			*(resi+i) = *(x+i)/(*(mu+i));

			//for omega:
			*(dmudtheta + i) = 1;
			for(k = 1; k <= *(order+1); k++) *(dmudtheta + i) += *(par + *order + k) * (*(dmudtheta + i - k));
			//for alpha:
			for(j = 1; j <= *order; j++){
				*(dmudtheta + j * (*N) + i) = *(x + i - 1);
				for(k = 1; k <= *(order+1); k++) *(dmudtheta + j * (*N) + i) +=	*(par + *order + k) * (*(dmudtheta + j*(*N) + i - k));
			}
			//for beta:
			for(j = *order + 1; j <= *(order+1) + *order; j++){
				*(dmudtheta + j*(*N) + i) = *(mu + i - 1);
				for(k = 1; k <= *(order+1); k++) *(dmudtheta + j*(*N) + i) += *(par + *order + k) * (*(dmudtheta + j*(*N) + i - k));
			}

		}
		startIndex = stopIndex;
		stopIndex = min(stopIndex + maxpq, *N);
	} while (stopIndex != *N);

}
//END---getdmudtheta_ACD----------------------------//
