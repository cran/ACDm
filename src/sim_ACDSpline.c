#include "header.h" 
#include <Rinternals.h>
#include <Rdefines.h>


//START---sim_ACDSpline----------------------------//
SEXP sim_ACDSpline(SEXP N,
		SEXP par,
		SEXP order,
		SEXP startX,
		SEXP startMu,
		SEXP e,
		SEXP Nburn,
		SEXP startTime,
		SEXP endTime,
		SEXP knots,
		SEXP constant,
		SEXP lin,
		SEXP sq,
		SEXP cube,
		SEXP splineNewDay){



	PROTECT(N = AS_INTEGER(N));
	PROTECT(par = AS_NUMERIC(par));
	PROTECT(order = AS_INTEGER(order));
	PROTECT(startX = AS_NUMERIC(startX));
	PROTECT(e = AS_NUMERIC(e));
	PROTECT(Nburn = AS_INTEGER(Nburn));
	PROTECT(startTime = AS_NUMERIC(startTime));
	PROTECT(endTime = AS_NUMERIC(endTime));
	PROTECT(splineNewDay = AS_INTEGER(splineNewDay));
	PROTECT(knots = AS_NUMERIC(knots));
	PROTECT(constant = AS_NUMERIC(constant));
	PROTECT(lin = AS_NUMERIC(lin));
	PROTECT(sq = AS_NUMERIC(sq));
	PROTECT(cube = AS_NUMERIC(cube));


	double *xpar, *xstartX, *xstartMu, *xe, *xdur, *xdurT, *xtime, timeTemp, *xstartTime, *xendTime;
	double *xknots, *xconstant, *xlin, *xsq, *xcube, adjFactor;
	int i = 0, j = 1, row = 0, stopCode = 0;
	int *xorder, *xsplineNewDay, *xday, *xN, *xNburn;
	int Nstart = LENGTH(startX), dayTemp = 0, knotPos = 0, lastPos = LENGTH(knots)-1;
	xpar = NUMERIC_POINTER(par); xorder = INTEGER_POINTER(order);
	xstartX = NUMERIC_POINTER(startX); xstartMu = NUMERIC_POINTER(startMu);
	xe = NUMERIC_POINTER(e); xstartTime = NUMERIC_POINTER(startTime); xendTime = NUMERIC_POINTER(endTime);
	xsplineNewDay = INTEGER_POINTER(splineNewDay);
	xknots = NUMERIC_POINTER(knots); xconstant = NUMERIC_POINTER(constant); xlin = NUMERIC_POINTER(lin); xsq = NUMERIC_POINTER(sq);
	xcube = NUMERIC_POINTER(cube); xN = INTEGER_POINTER(N); xNburn = INTEGER_POINTER(Nburn);

	SEXP dur, durT, timeC, day, list;
	double mu[*xN + *xNburn];
	double xTemp[*xN + *xNburn];
	PROTECT(dur = NEW_NUMERIC(*xN));
	xdur = NUMERIC_POINTER(dur) ;
	PROTECT(durT = NEW_NUMERIC(*xN));
	xdurT = NUMERIC_POINTER(durT) ;
	PROTECT(timeC = NEW_NUMERIC(*xN));
	xtime = NUMERIC_POINTER(timeC);
	PROTECT(day = NEW_INTEGER(*xN));
	xday = INTEGER_POINTER(day);
	PROTECT(list = NEW_LIST(4));

	do {
		for (i = 0; i < Nstart; i++){ //for the start values
			xTemp[i] = xstartX[i];
			mu[i] = xstartMu[i];
		}
		for (; i < *xNburn; i++){ //for the burn
			mu[i] = xpar[0]; //adds the constant
			for(j = 1; j <= xorder[0]; j++) mu[i] += xpar[j]*xTemp[i-j]; //adds the lagged durations
			for(j = 1; j <= xorder[1]; j++) mu[i] += xpar[j+xorder[0]]*mu[i-j]; //adds the lagged mus
			xTemp[i] = mu[i]*xe[i];
		}
		mu[i] = xpar[0]; //adds the constant
		for(j = 1; j <= xorder[0]; j++) mu[i] += xpar[j]*xTemp[i-j]; //adds the lagged durations
		for(j = 1; j <= xorder[1]; j++) mu[i] += xpar[j+xorder[0]]*mu[i-j]; //adds the lagged mus
		xTemp[i] = mu[i]*xe[i];
		xdur[row] = xTemp[i];

		knotPos = xsplineNewDay[dayTemp%5];
		timeTemp = xstartTime[0] - xknots[xsplineNewDay[dayTemp%5]];
		adjFactor = xconstant[xsplineNewDay[dayTemp%5]] +
				timeTemp * xlin[xsplineNewDay[dayTemp%5]] +
				pow(timeTemp, 2) * xsq[xsplineNewDay[dayTemp%5]] +
				pow(timeTemp, 3) * xcube[xsplineNewDay[dayTemp%5]];

		xdurT[row] = xTemp[i] * adjFactor;
		xtime[row] = xstartTime[0] + xdurT[row];
		xday[row] = dayTemp;
		i++; row++;
		do{ //for the real series
			stopCode = 0;
			mu[i] = xpar[0]; //adds the constant
			for(j = 1; j <= xorder[0]; j++) mu[i] += xpar[j]*xTemp[i-j]; //adds the lagged durations
			for(j = 1; j <= xorder[1]; j++) mu[i] += xpar[j+xorder[0]]*mu[i-j]; //adds the lagged mus
			xTemp[i] = mu[i]*xe[i];
			xdur[row] = xTemp[i];

			if(dayTemp%5 != 4){
				while((knotPos < xsplineNewDay[(dayTemp+1)%5] - 1) &&
						(xtime[row-1] > xknots[knotPos+1])) knotPos++;
				timeTemp = xtime[row-1] - xknots[knotPos];
				if(xtime[row-1] < xknots[xsplineNewDay[dayTemp%5]]){ //outside the knots the function is linear
					adjFactor = xconstant[knotPos] +
							timeTemp * xlin[knotPos];
				} else{
					adjFactor = xconstant[knotPos] +
							timeTemp * xlin[knotPos] +
							pow(timeTemp, 2) * xsq[knotPos] +
							pow(timeTemp, 3) * xcube[knotPos];
				}
			} else{
				while((knotPos < lastPos) &&
						(xtime[row-1] > xknots[knotPos+1])) knotPos++;
				timeTemp = xtime[row-1] - xknots[knotPos];
				if(xtime[row-1] < xknots[xsplineNewDay[dayTemp%5]]){ //outside the knots the function is linear
					adjFactor = xconstant[knotPos] +
							timeTemp * xlin[knotPos];
				} else{
					adjFactor = xconstant[knotPos] +
							timeTemp * xlin[knotPos] +
							pow(timeTemp, 2) * xsq[knotPos] +
							pow(timeTemp, 3) * xcube[knotPos];
				}
			}


			xdurT[row] = xTemp[i] * adjFactor;
			xtime[row] = xtime[row-1] + xdurT[row];
			xday[row] = dayTemp;

			if(xtime[row] > xendTime[0]){
				stopCode = 2;
				row--;
			}
			else if(row >= *xN-1) stopCode = 1;
			row++;
			i++;
		} while(stopCode == 0);
		dayTemp++;
	} while (stopCode != 1);

	SET_VECTOR_ELT(list, 0, day);
	SET_VECTOR_ELT(list, 1, timeC);
	SET_VECTOR_ELT(list, 2, durT);
	SET_VECTOR_ELT(list, 3, dur);

	UNPROTECT(19);
	return(list);
}
//END---sim_ACDSpline----------------------------//
