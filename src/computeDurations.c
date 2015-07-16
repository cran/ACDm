#include "header.h" 

//START---computeDurationsSubSec----------------------------//
void computeDurationsSubSec(int *y,
		int *M,
		int *d,
		int *h,
		int *m,
		double *s,
		int *yDur,
		int *MDur,
		int *dDur,
		int *hDur,
		int *mDur,
		double *sDur,
		int *vol,
		double *price,
		int *volDur,
		double *priceDur,
		int *Ntrans,
		int *dur,
		int *Ntime,
		int *Ndur,
		double *open,
		double *close,
		int *durType,
		int *zeroDurHandeling,
		double *priceChange,
		int *culmVol){

	int i, j = 0, lastMonth = 0, lastDay = 0, lastYear = 0, tempNtrans = 1,  tempDur = 0, tempVol = 0;
	double tempTurnover = 0, lastPrice = 0, secOfDay = 0, lastSecOfDay = 0;

	if((*durType == 1) && (*zeroDurHandeling == 1)){ // trade durations with 0s removed/aggregated
		for (i = 0; i < *Ntime-1; i++) {
			if(!(lastYear == *(y+i) && lastMonth == *(M+i) && lastDay == *(d+i))){
				while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) <= *open)	i++;
				lastYear = *(y+i-1); lastMonth = *(M+i-1); lastDay = *(d+i-1);
				lastSecOfDay = *open;
				tempVol = 0; tempTurnover = 0; tempNtrans = 1;
			}
			secOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
			if(secOfDay <= *close){
				if(secOfDay == *(h+i+1)*3600 + *(m+i+1)*60 + *(s+i+1)){ //next duration is zero
					tempVol += *(vol+i);
					tempTurnover += *(vol+i) * (*(price+i));
					tempNtrans++;
					tempDur = secOfDay - lastSecOfDay;
				}
				else if(((*(h+i-1)*3600 + *(m+i-1)*60 + *(s+i-1))) == secOfDay){ //this duration is zero but not next
					tempVol += (*(vol+i));
					tempTurnover += (*(vol+i)) * (*(price+i));
					//adds a row to the new duration object:
					*(yDur+j) = *(y + i);
					*(MDur+j) = *(M + i);
					*(dDur+j) = *(d + i);
					*(hDur+j) = *(h + i);
					*(mDur+j) = *(m + i);
					*(sDur+j) = *(s + i);
					*(volDur+j) = tempVol;
					*(priceDur+j) = (tempTurnover/tempVol);
					*(dur+j) = tempDur;
					*(Ntrans+j) = tempNtrans;

					lastSecOfDay = secOfDay;
					j++;
					tempVol = 0;
					tempTurnover = 0;
					tempNtrans = 1;
				}
				else{ //neither the current duration or the next is zero
					*(yDur+j) = *(y + i);
					*(MDur+j) = *(M + i);
					*(dDur+j) = *(d + i);
					*(hDur+j) = *(h + i);
					*(mDur+j) = *(m + i);
					*(sDur+j) = *(s + i);
					*(volDur+j) = *(vol + i);
					*(priceDur+j) = *(price + i);
					*(dur+j) = secOfDay - lastSecOfDay;
					*(Ntrans+j) = 1;

					lastSecOfDay = secOfDay;
					j++;
				}
			}
		}

		//for the last observation:
		secOfDay = *(h+*Ntime) * 3600 + *(m+*Ntime) * 60 + *(s+*Ntime);
		if(tempNtrans > 1){ //zero duration
			tempVol += *(vol+*Ntime);
			tempTurnover += *(vol+*Ntime)*(*(price+*Ntime));
			//adds a row to the new duration object:
			*(yDur+j) = *(y + *Ntime);
			*(MDur+j) = *(M + *Ntime);
			*(dDur+j) = *(d + *Ntime);
			*(hDur+j) = *(h + *Ntime);
			*(mDur+j) = *(m + *Ntime);
			*(sDur+j) = *(s + *Ntime);
			*(volDur+j) = tempVol;
			*(priceDur+j) = tempTurnover/tempVol;
			*(dur+j) = secOfDay - (*(h+*Ntime-tempNtrans)*3600 + *(m+*Ntime-tempNtrans)*60 + *(s+*Ntime-tempNtrans));
			*(Ntrans+j) = tempNtrans;
		}
		else { //not zero duration
			*(yDur+j) = *(y + *Ntime);
			*(MDur+j) = *(M + *Ntime);
			*(dDur+j) = *(d + *Ntime);
			*(hDur+j) = *(h + *Ntime);
			*(mDur+j) = *(m + *Ntime);
			*(sDur+j) = *(s + *Ntime);
			*(volDur+j) = *(vol + *Ntime);
			*(priceDur+j) = *(price + *Ntime);
			*(dur+j) = secOfDay - (*(h + *Ntime - 1) * 3600 + *(m + *Ntime - 1) * 60 + *(s + *Ntime - 1));
			*(Ntrans+j) = 1;
		}
	}
	else if((*durType == 1) && (*zeroDurHandeling == 0)){ // trade durations with 0s kept
		for (i = 0; i < *Ntime; i++) {
							if(!(lastYear == *(y+i) && lastMonth == *(M+i) && lastDay == *(d+i))){ //if new day
								while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) <= *open)	i++; //removes transactions before opening
								lastYear = *(y+i-1); lastMonth = *(M+i-1); lastDay = *(d+i-1);
								lastSecOfDay = *open;
								tempNtrans = 1;
							}
							secOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
							if(secOfDay <= *close){
								//adds a row to the new duration object:
								*(yDur+j) = *(y + i);
								*(MDur+j) = *(M + i);
								*(dDur+j) = *(d + i);
								*(hDur+j) = *(h + i);
								*(mDur+j) = *(m + i);
								*(sDur+j) = *(s + i);
								*(volDur+j) = *(vol + i);
								*(priceDur+j) = *(price + i);
								*(dur+j) = secOfDay - lastSecOfDay;
								lastSecOfDay = secOfDay;
								j++;
							}
		}
	}
	else if(*durType == 2){ // price durations
		if(*zeroDurHandeling == 1){ //0s should  be aggregated
			for (i = 0; i < *Ntime; i++) {
				if(!(lastYear == *(y+i) && lastMonth == *(M+i) && lastDay == *(d+i))){
					while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) < *open)	i++;
					if((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open){
						while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open)	i++;
						lastPrice = *(price+i-1); //if there were transaction the first second of the trade day, the last price during that second will be used
						lastSecOfDay = *open;
					}
					else{
						lastPrice = *(price+i); //if there were no transactions the first second of the day, the first transaction of the day will be used
						lastSecOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
						i++;
					}
					lastYear = *(y+i-1); lastMonth = *(M+i-1); lastDay = *(d+i-1);
					tempVol = 0; tempTurnover = 0; tempNtrans = 0;
				}
				secOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);

				if(secOfDay <= *close){
					tempVol += *(vol+i);
					tempNtrans++;

					if(fabs(*(price+i) - lastPrice) >= *priceChange && secOfDay != lastSecOfDay){ //price difference is big enough and duration is not zero
						//adds a row to the new duration object:
						*(yDur+j) = *(y + i);
						*(MDur+j) = *(M + i);
						*(dDur+j) = *(d + i);
						*(hDur+j) = *(h + i);
						*(mDur+j) = *(m + i);
						*(sDur+j) = *(s + i);
						*(volDur+j) = tempVol;
						*(priceDur+j) = *(price+i);
						*(dur+j) = secOfDay - lastSecOfDay;
						*(Ntrans+j) = tempNtrans;

						lastSecOfDay = secOfDay;
						lastPrice = *(price+i);
						j++;
						tempVol = 0;
						tempTurnover = 0;
						tempNtrans = 0;
					}
				}
			}
		}
		if(*zeroDurHandeling == 0){ //0s should not be aggregated
					for (i = 0; i < *Ntime; i++) {
						if(!(lastYear == *(y+i) && lastMonth == *(M+i) && lastDay == *(d+i))){
							while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) < *open)	i++;
							if((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open){
								while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open)	i++;
								lastPrice = *(price+i-1); //if there were transaction the first second of the trade day, the last price during that second will be used
								lastSecOfDay = *open;
							}
							else{
								lastPrice = *(price+i); //if there were no transactions the first second of the day, the first transaction of the day will be used
								lastSecOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
								i++;
							}
							lastYear = *(y+i-1); lastMonth = *(M+i-1); lastDay = *(d+i-1);
							tempVol = 0; tempTurnover = 0; tempNtrans = 0;
						}
						secOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);

						if(secOfDay <= *close){
							tempVol += *(vol+i);
							tempNtrans++;

							if(fabs(*(price+i) - lastPrice) >= *priceChange){ //price difference is big enough
								//adds a row to the new duration object:
								*(yDur+j) = *(y + i);
								*(MDur+j) = *(M + i);
								*(dDur+j) = *(d + i);
								*(hDur+j) = *(h + i);
								*(mDur+j) = *(m + i);
								*(sDur+j) = *(s + i);
								*(volDur+j) = tempVol;
								*(priceDur+j) = *(price+i);
								*(dur+j) = secOfDay - lastSecOfDay;
								*(Ntrans+j) = tempNtrans;

								lastSecOfDay = secOfDay;
								lastPrice = *(price+i);
								j++;
								tempVol = 0;
								tempTurnover = 0;
								tempNtrans = 0;
							}
						}
					}
				}
	}
	else if(*durType == 3){ // volume durations
			if(*zeroDurHandeling == 1){ //0s should  be aggregated
				for (i = 0; i < *Ntime; i++) {
					if(!(lastYear == *(y+i) && lastMonth == *(M+i) && lastDay == *(d+i))){
						while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) < *open)	i++;
						if((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open){
							while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open)	i++;
							//if there were transaction the first second of the trade day, the last price during that second will be used
							lastSecOfDay = *open;
						}
						else{
							//if there were no transactions the first second of the day, the first transaction of the day will be used
							lastSecOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
							i++;
						}
						lastYear = *(y+i-1); lastMonth = *(M+i-1); lastDay = *(d+i-1);
						tempVol = 0; tempTurnover = 0; tempNtrans = 0;
					}
					secOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);

					if(secOfDay <= *close){
						tempVol += *(vol+i);
						tempNtrans++;

						if(tempVol >= *culmVol && secOfDay != lastSecOfDay){ //Cumulated volume is large enough and duration is not zero
							//adds a row to the new duration object:
							*(yDur+j) = *(y + i);
							*(MDur+j) = *(M + i);
							*(dDur+j) = *(d + i);
							*(hDur+j) = *(h + i);
							*(mDur+j) = *(m + i);
							*(sDur+j) = *(s + i);
							*(volDur+j) = tempVol;
							*(priceDur+j) = *(price+i);
							*(dur+j) = secOfDay - lastSecOfDay;
							*(Ntrans+j) = tempNtrans;

							lastSecOfDay = secOfDay;
							j++;
							tempVol = 0;
							tempTurnover = 0;
							tempNtrans = 0;
						}
					}
				}
			} else if(*zeroDurHandeling == 0){ //0s should  not be aggregated
				for (i = 0; i < *Ntime; i++) {
					if(!(lastYear == *(y+i) && lastMonth == *(M+i) && lastDay == *(d+i))){
						while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) < *open)	i++;
						if((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open){
							while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open)	i++;
							//if there were transaction the first second of the trade day, the last price during that second will be used
							lastSecOfDay = *open;
						}
						else{
							//if there were no transactions the first second of the day, the first transaction of the day will be used
							lastSecOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
							i++;
						}
						lastYear = *(y+i-1); lastMonth = *(M+i-1); lastDay = *(d+i-1);
						tempVol = 0; tempTurnover = 0; tempNtrans = 0;
					}
					secOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);

					if(secOfDay <= *close){
						tempVol += *(vol+i);
						tempNtrans++;

						if(tempVol >= *culmVol){ //Cumulated volume is large enough
							//adds a row to the new duration object:
							*(yDur+j) = *(y + i);
							*(MDur+j) = *(M + i);
							*(dDur+j) = *(d + i);
							*(hDur+j) = *(h + i);
							*(mDur+j) = *(m + i);
							*(sDur+j) = *(s + i);
							*(volDur+j) = tempVol;
							*(priceDur+j) = *(price+i);
							*(dur+j) = secOfDay - lastSecOfDay;
							*(Ntrans+j) = tempNtrans;

							lastSecOfDay = secOfDay;
							j++;
							tempVol = 0;
							tempTurnover = 0;
							tempNtrans = 0;
						}
					}
				}
			}
		}

	*Ndur = j;
}
//END---computeDurationsSubSec----------------------------//

//START---computeDurationsShort----------------------------//
void computeDurationsShort(int *y,
		int *M,
		int *d,
		int *h,
		int *m, //5
		int *s,
		int *yDur,
		int *MDur,
		int *dDur,
		int *hDur, //10
		int *mDur,
		int *sDur,
		int *dur,
		int *Ndur,
		int *Ntrans, //15
		int *Ntime,
		int *open,
		int *close,
		int *zeroDurHandeling){ //19

	int i, j = 0, lastMonth = 0, lastDay = 0, lastYear = 0, tempNtrans = 1, sekOfDay = 0, lastSekOfDay = 0, tempDur = 0;


	if(*zeroDurHandeling == 1){ //0s removed/aggregated
		for (i = 0; i < *Ntime-1; i++) {
			if(!(lastYear == *(y+i) && lastMonth == *(M+i) && lastDay == *(d+i))){ //if new day
				while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) <= *open)	i++; //removes transactions before opening
				lastYear = *(y+i-1); lastMonth = *(M+i-1); lastDay = *(d+i-1);
				lastSekOfDay = *open;
				tempNtrans = 1;
			}
			sekOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
			if(sekOfDay <= *close){
				if(sekOfDay == *(h+i+1)*3600 + *(m+i+1)*60 + *(s+i+1)){ //next duration is zero
					tempNtrans++;
					tempDur = sekOfDay - lastSekOfDay;
				}
				else if(((*(h+i-1)*3600 + *(m+i-1)*60 + *(s+i-1))) == sekOfDay){ //this duration is zero but not next
					//adds a row to the new duration object:
					*(yDur+j) = *(y + i);
					*(MDur+j) = *(M + i);
					*(dDur+j) = *(d + i);
					*(hDur+j) = *(h + i);
					*(mDur+j) = *(m + i);
					*(sDur+j) = *(s + i);
					*(dur+j) = tempDur;
					*(Ntrans+j) = tempNtrans;

					lastSekOfDay = sekOfDay;
					j++;
					tempNtrans = 1;
				}
				else{ //neither the current duration nor the next is zero
					*(yDur+j) = *(y + i);
					*(MDur+j) = *(M + i);
					*(dDur+j) = *(d + i);
					*(hDur+j) = *(h + i);
					*(mDur+j) = *(m + i);
					*(sDur+j) = *(s + i);
					*(dur+j) = sekOfDay - lastSekOfDay;
					*(Ntrans+j) = 1;

					lastSekOfDay = sekOfDay;
					j++;
				}
			}
		}

		//for the last observation:
		sekOfDay = *(h+*Ntime) * 3600 + *(m+*Ntime) * 60 + *(s+*Ntime);
		if(tempNtrans > 1){ //zero duration
			//adds a row to the new duration object:
			*(yDur+j) = *(y + *Ntime);
			*(MDur+j) = *(M + *Ntime);
			*(dDur+j) = *(d + *Ntime);
			*(hDur+j) = *(h + *Ntime);
			*(mDur+j) = *(m + *Ntime);
			*(sDur+j) = *(s + *Ntime);
			*(dur+j) = sekOfDay - (*(h+*Ntime-tempNtrans)*3600 + *(m+*Ntime-tempNtrans)*60 + *(s+*Ntime-tempNtrans));
			*(Ntrans+j) = tempNtrans;
		} else { //not zero duration
			*(yDur+j) = *(y + *Ntime);
			*(MDur+j) = *(M + *Ntime);
			*(dDur+j) = *(d + *Ntime);
			*(hDur+j) = *(h + *Ntime);
			*(mDur+j) = *(m + *Ntime);
			*(sDur+j) = *(s + *Ntime);
			*(dur+j) = sekOfDay - (*(h + *Ntime - 1) * 3600 + *(m + *Ntime - 1) * 60 + *(s + *Ntime - 1));
			*(Ntrans+j) = 1;
		}
	} else if(*zeroDurHandeling == 0){ //zeroes not removed
		for (i = 0; i < *Ntime; i++) {
					if(!(lastYear == *(y+i) && lastMonth == *(M+i) && lastDay == *(d+i))){ //if new day
						while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) <= *open)	i++; //removes transactions before opening
						lastYear = *(y+i-1); lastMonth = *(M+i-1); lastDay = *(d+i-1);
						lastSekOfDay = *open;
						tempNtrans = 1;
					}
					sekOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
					if(sekOfDay <= *close){
						//adds a row to the new duration object:
						*(yDur+j) = *(y + i);
						*(MDur+j) = *(M + i);
						*(dDur+j) = *(d + i);
						*(hDur+j) = *(h + i);
						*(mDur+j) = *(m + i);
						*(sDur+j) = *(s + i);
						*(dur+j) = sekOfDay - lastSekOfDay;
						lastSekOfDay = sekOfDay;
						j++;
					}
		}
	}

	*Ndur = j;
}
//END---computeDurationsShort----------------------------//
