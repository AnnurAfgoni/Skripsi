#include<stdio.h>
#include<math.h>

#define M_PI 3.14159265358979323846

int main(){
	
	int i, j;
	int m, n, p, q;
	int n_bar = 10, n_col = 10;
	
	int NX, NY, LX, LY;
	NX = 100; NY = 100; LX = 10; LY = 10;
	
	float dx, dy;
	dx = (float)LX/(float)NX;
	dy = (float)LY/(float)NY;
	
	float x0, y0;
	x0 = (float)LX/2.0;
	y0 = (float)LY/2.0;
	
	double element;
	
	double x[NX+1][NY+1], y[NX+1][NY+1];
	double V[NX+1][NY+1];
	double psi_kiri[NX+1][NY+1], psi_kanan[NX+1][NY+1];
	double d2psix[NX+1][NY+1], d2psiy[NX+1][NY+1];
	
	double represent[n_bar*n_col][n_bar*n_col];
	
	// Data input
	for(i=0; i<NX+1; i++){
	for(j=0; j<NY+1; j++){
		x[i][j] = (float)j*dx;
		y[i][j] = (float)i*dy;
	}}
	/*
	for(i=0; i<NX+1; i++){
	for(j=0; j<NY+1; j++){
		printf("%.1lf ", y[i][j]);
	}printf("\n");}
	*/
	
	
	//Potensial OH
	for(i=0; i<NX+1; i++){
	for(j=0; j<NY+1; j++){
		V[i][j] = 0.5*(x[i][j] - x0)*(x[i][j] - x0) + 0.5*(y[i][j] - y0)*(y[i][j] - y0);
	}}
	/*
	for(i=0; i<NX+1; i++){
	for(j=0; j<NY+1; j++){
		printf("%.1lf ", V[i][j]);
	}printf("\n");}
	*/
	
	//Membuat Matriks
	for(m=1; m<=n_bar; m++){
	for(n=1; n<=n_col; n++){
	for(p=1; p<=n_bar; p++){
	for(q=1; q<=n_col; q++){
		
		element = 0.0;
		
		for(i=0; i<NX+1; i++){
		for(j=0; j<NY+1; j++){
			psi_kiri[i][j] = 2/sqrt((float)LX*(float)LY) * sin((float)m*M_PI*x[i][j]/(float)LX) * sin((float)n*M_PI*y[i][j]/(float)LX);
			psi_kanan[i][j] = 2/sqrt((float)LX*(float)LY) * sin((float)p*M_PI*x[i][j]/(float)LY) * sin((float)q*M_PI*y[i][j]/(float)LY);
		}}
		
		for(i=0; i<NX+1; i++){
		for(j=0; j<NY+1; j++){
			d2psix[i][j] = -0.5 * -((float)p*M_PI/(float)LX)*((float)p*M_PI/(float)LX) * 2/sqrt((float)LX*(float)LY) * sin((float)p*M_PI*x[i][j]/(float)LX) * sin((float)q*M_PI*x[i][j]/(float)LX);
			d2psiy[i][j] = -0.5 * -((float)q*M_PI/(float)LX)*((float)q*M_PI/(float)LX) * 2/sqrt((float)LX*(float)LY) * sin((float)p*M_PI*x[i][j]/(float)LY) * sin((float)q*M_PI*y[i][j]/(float)LY);
		}}
		
		for(i=1; i<NX; i++){
		for(j=1; j<NY; j++){
			element += psi_kiri[i][j] * (d2psix[i][j] + d2psiy[i][j]) + psi_kiri[i][j]*V[i][j]*psi_kanan[i][j];
		}}
		
		element *= dx*dy;
		
		printf("%d%d_%d%d = %lf \n", m, n, p, q, element);
		
	}}}}
	
	
	return 0;
}