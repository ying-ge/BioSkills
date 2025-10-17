#include<stdio.h>
#include<math.h>
#include<R.h>
#include<Rmath.h>
#include<Rinternals.h>


int csign(double X_i, int Censor_i, double X_j, int Censor_j) {

if (X_i > X_j) {
if ( (Censor_i==1) && (Censor_j==1) ) return(1);
if ( (Censor_i==1) && (Censor_j==0) ) return(0);
if ( (Censor_i==0) && (Censor_j==1) ) return(1);
if ( (Censor_i==0) && (Censor_j==0) ) return(0); }

if (X_i < X_j) {
if ( (Censor_i==1) && (Censor_j==1) ) return(-1);
if ( (Censor_i==1) && (Censor_j==0) ) return(-1);
if ( (Censor_i==0) && (Censor_j==1) ) return(0);
if ( (Censor_i==0) && (Censor_j==0) ) return(0); }

if (X_i == X_j) {
if ( (Censor_i==1) && (Censor_j==1) ) return(0);
if ( (Censor_i==1) && (Censor_j==0) ) return(-1);
if ( (Censor_i==0) && (Censor_j==1) ) return(1);
if ( (Censor_i==0) && (Censor_j==0) ) return(0); }

return 0; // -1 should cause the application to throw an exception and crash, help track down errors
}



void TauXX(double *timeX, int *statusX, int *nn, double *output) {	
								// timeX is a continuous variable
								// statusX 1 indicates event, 0 indicates censored
int i, j, nobs = *nn;
double est_tXX = 0;

for (i=0; i<nobs; i++) {
for (j=0; j<nobs; j++) {  if (j==i) continue;
est_tXX += csign(timeX[i],statusX[i],timeX[j],statusX[j])*csign(timeX[i],statusX[i],timeX[j],statusX[j]);
}}

*output = (double) est_tXX/nobs/(nobs-1);
}



void TauXY(double *timeX, int *statusX, double *scoreY, int *nn, double *output) {

int i, j, nobs = *nn;
double est_tXY = 0;

for (i=0; i<nobs; i++) {
for (j=0; j<nobs; j++) {  if (j==i) continue;
est_tXY += csign(timeX[i],statusX[i],timeX[j],statusX[j])*sign(scoreY[i]-scoreY[j]);
}}

*output = (double) est_tXY/nobs/(nobs-1);
}


void VarTauXX(double *timeX, int *statusX, int *nn, double *output) {

int i, j, nobs = *nn;
double temp = 0;
double temp_s1 = 0, temp_s2 = 0, temp_s3 = 0;

for (i=0; i<nobs; i++) {  double temp_s1_j = 0, temp_s3_j = 0;
for (j=0; j<nobs; j++) {  if (j==i) continue;
temp = csign(timeX[i],statusX[i],timeX[j],statusX[j])*csign(timeX[i],statusX[i],timeX[j],statusX[j]);
temp_s1_j += temp;
temp_s3_j += temp*temp;
}

temp_s1 += 4*temp_s1_j*temp_s1_j;
temp_s3 += 2*temp_s3_j;
temp_s2 += temp_s1_j;
}

temp_s2 = (double) temp_s2/nobs/(nobs-1)*temp_s2*(2*nobs-3)*(-2);

*output = (double) (temp_s1-temp_s3+temp_s2)/nobs/(nobs-1)/(nobs-2)/(nobs-3);
}



void VarTauXY(double *timeX, int *statusX, double *scoreY, int *nn, double *output) {

int i, j, nobs = *nn;
double temp = 0;
double temp_s1 = 0, temp_s2 = 0, temp_s3 = 0;

for (i=0; i<nobs; i++) {  double temp_s1_j = 0, temp_s3_j = 0;
for (j=0; j<nobs; j++) {  if (j==i) continue;
temp = csign(timeX[i],statusX[i],timeX[j],statusX[j])*sign(scoreY[i]-scoreY[j]);
temp_s1_j += temp;
temp_s3_j += temp*temp;
}

temp_s1 += 4*temp_s1_j*temp_s1_j;
temp_s3 += 2*temp_s3_j;
temp_s2 += temp_s1_j;
}

temp_s2 = (double) temp_s2/nobs/(nobs-1)*temp_s2*(2*nobs-3)*(-2);

*output = (double) (temp_s1-temp_s3+temp_s2)/nobs/(nobs-1)/(nobs-2)/(nobs-3);
}





void CovTauXXXY(double *timeX, int *statusX, double *scoreY, int *nn, double *output) {

int i, j, nobs = *nn;
double temp_XX = 0, temp_XY = 0;
double temp_s1 = 0, temp_s2 = 0, temp_s3 = 0, temp_s4 = 0;

for (i=0; i<nobs; i++) {  double temp_sXX_j = 0, temp_sXY_j = 0, temp_sXXXY_j = 0;
for (j=0; j<nobs; j++) {  if (j==i) continue;
temp_XX = csign(timeX[i],statusX[i],timeX[j],statusX[j])*csign(timeX[i],statusX[i],timeX[j],statusX[j]);
temp_XY = csign(timeX[i],statusX[i],timeX[j],statusX[j])*sign(scoreY[i]-scoreY[j]);
temp_sXX_j   += temp_XX;
temp_sXY_j   += temp_XY;
temp_sXXXY_j += temp_XX*temp_XY;
}

temp_s1 += 4*temp_sXX_j*temp_sXY_j;
temp_s3 += 2*temp_sXXXY_j;

temp_s2 += temp_sXX_j;
temp_s4 += temp_sXY_j;
}

*output = (double) (temp_s1-temp_s3+(2*nobs-3)*(-2)*temp_s2*temp_s4/nobs/(nobs-1))/nobs/(nobs-1)/(nobs-2)/(nobs-3);
}




void CovTauXYXZ(double *timeX, int *statusX, double *scoreY, double *scoreZ, int *nn, double *output) {

int i, j, nobs = *nn;
double temp_XY = 0, temp_XZ = 0;
double temp_s1 = 0, temp_s2 = 0, temp_s3 = 0, temp_s4 = 0;

for (i=0; i<nobs; i++) {  double temp_sXY_j = 0, temp_sXZ_j = 0, temp_sXYXZ_j = 0;
for (j=0; j<nobs; j++) {  if (j==i) continue;
temp_XY = csign(timeX[i],statusX[i],timeX[j],statusX[j])*sign(scoreY[i]-scoreY[j]);
temp_XZ = csign(timeX[i],statusX[i],timeX[j],statusX[j])*sign(scoreZ[i]-scoreZ[j]);
temp_sXY_j   += temp_XY;
temp_sXZ_j   += temp_XZ;
temp_sXYXZ_j += temp_XY*temp_XZ;
}

temp_s1 += 4*temp_sXY_j*temp_sXZ_j;
temp_s3 += 2*temp_sXYXZ_j;

temp_s2 += temp_sXY_j;
temp_s4 += temp_sXZ_j;
}

*output = (double) (temp_s1-temp_s3+(2*nobs-3)*(-2)*temp_s2*temp_s4/nobs/(nobs-1))/nobs/(nobs-1)/(nobs-2)/(nobs-3);
}



