#include<stdio.h>
#include<math.h>

//#define M_PI 3.141592653589793

int main()
{
	int i, j, m, n, p, q, ii, jj;
	int NX = 100, NY = 100, LX = 10, LY = 10;
	
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
	
	// Data input
	for(i=0; i<NX+1; i++){
	for(j=0; j<NY+1; j++){
		x[i][j] = (float)j*dx;
		y[i][j] = (float)i*dy;
	}}
	
	//Potensial OH
	for(i=0; i<NX+1; i++){
	for(j=0; j<NY+1; j++){
		V[i][j] = 0.5*(x[i][j] - x0)*(x[i][j] - x0) + 0.5*(y[i][j] - y0)*(y[i][j] - y0);
	}}
	
	// indexing
	i = 0;
	int ind[100][2];
	for(m = 1; m <= 10; m++){
	for(n = 1; n <= 10; n++){
		ind[i][0] = m;
		ind[i][1] = n;
		i++;
	}}	
	
	// konstanta
	double pix, piy;
	double a = 2/sqrt((float)LX*(float)LY);
	pix = M_PI/(double)LX;
	piy = M_PI/(double)LY;
	
	double rep_matriks[NX][NY];
	
	for(ii=0; ii<NX; ii++){
		m = ind[ii][0];
		n = ind[ii][1];
		
		for(i=0; i<NX+1; i++){
		for(j=0; j<NY+1; j++){
			psi_kiri[i][j] = a * sin((float)m*pix*x[i][j]) * sin((float)n*piy*y[i][j]);
		}}
		
		for(jj=0; jj<NY; jj++){					
			p = ind[jj][0];
			q = ind[jj][1];			
			
			for(i=0; i<NX+1; i++){
			for(j=0; j<NY+1; j++){
				psi_kanan[i][j] = a * sin((float)p*pix*x[i][j]) * sin((float)q*piy*y[i][j]);
			}}
			
			for(i=0; i<NX+1; i++){
			for(j=0; j<NY+1; j++){
				d2psix[i][j] = -0.5 * -((float)p*piy)*((float)p*piy) * a * sin((float)p*piy*x[i][j]) * sin((float)q*piy*y[i][j]);
				d2psiy[i][j] = -0.5 * -((float)p*piy)*((float)p*piy) * a * sin((float)p*piy*x[i][j]) * sin((float)q*piy*y[i][j]);
			}}
			
			element = 0.0;
			for(i=1; i<NX; i++){
			for(j=1; j<NY; j++){
				element += psi_kiri[i][j]*(d2psix[i][j]+d2psiy[i][j]) + psi_kiri[i][j]*V[i][j]*psi_kanan[i][j];
			}}
			
			element *= dx*dy;
			
			// printf("%d%d %d%d = %.15lf\n",m,n,p,q,element);
			
			rep_matriks[ii][jj] = element;
		}
	}
	
	FILE *outfile;
	
	// simpan
	outfile = fopen("matriks_c.csv", "w");
	for(i=0; i<NX; i++){
		for(j=0; j<NY; j++){
			fprintf(outfile, "%.20lf,", rep_matriks[i][j]);
		}
		fprintf(outfile, "\n");
	}
	
}