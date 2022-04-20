#include<cmath>
#include<iostream>
#include<Eigen/Dense>

using namespace std;
using namespace Eigen;

int main(){
	
	int i, j, m, n, p, q;
	int n_bar = 10, n_col = 10;
	
	int NX = 100, NY = 100, LX = 10, LY = 10;
	
	float dx = float(LX)/float(NX), dy = float(LY)/float(NY);
	
	float x0, y0;
	x0 = float(LX)/2.0;
	y0 = float(LY)/2.0;
	MatrixXd xo, yo;
	xo = MatrixXd::Constant(NX+1, NY+1, x0);
	yo = MatrixXd::Constant(NX+1, NY+1, y0);

	double element;
	
	// data input
	MatrixXd x(NX+1, NY+1), y(NX+1, NY+1);
	VectorXd isi = VectorXd::LinSpaced(100+1, 0, 10);
	for(i=0; i<NX+1; i++){
		x.row(i) = isi;
		y.col(i) = isi;
	}
	
	// Potensial OH
	MatrixXd V(NX+1, NY+1);
	for(i=0; i<NX+1; i++){
	for(j=0; j<NY+1; j++){
		V(i,j) = 0.5*(x(i,j) - xo(i,j))*(x(i,j) - xo(i,j)) + 0.5*(y(i,j) - yo(i,j))*(y(i,j) - yo(i,j));	
	}}	
	
	//Membuat Matriks
	MatrixXd psi_kiri(NX+1, NY+1), psi_kanan(NX+1, NY+1);
	MatrixXd d2psix(NX+1, NY+1), d2psiy(NX+1, NY+1);
	MatrixXd rep_matriks(NX+1, NY+1);
	
	double a = 2/sqrt(LX*LY);
	double psix = M_PI/LX, psiy = M_PI/LY;
	
	for(m=1; m<=n_bar; m++){
	for(n=1; n<=n_col; n++){
		
		for(i=0; i<NX+1; i++){
		for(j=0; j<NY+1; j++){
			psi_kiri(i, j) = a * sin(m*psix*x(i,j)) * sin(n*psiy*y(i,j));
		}}
		
		for(p=1; p<=n_bar; p++){
		for(q=1; q<=n_col; q++){
			
			element = 0.0;
			
			for(i=0; i<NX+1; i++){
			for(j=0; j<NY+1; j++){
				psi_kanan(i,j) = a * sin(p*psiy*x(i,j)) * sin(q*psiy*y(i,j));
			}}
			
			for(i=0; i<NX+1; i++){
			for(j=0; j<NY+1; j++){
				d2psix(i,j) = -0.5 * -(p*psiy)*(p*psiy) * a * sin(p*psiy*x(i,j)) * sin(q*psiy*y(i,j));
				d2psiy(i,j) = -0.5 * -(q*psiy)*(q*psiy) * a * sin(p*psiy*x(i,j)) * sin(q*psiy*y(i,j));
			}}
			
			for(i=1; i<NX; i++){
			for(j=1; j<NY; j++){
				element += psi_kiri(i,j)*(d2psix(i,j)+d2psiy(i,j)) + psi_kiri(i,j)*V(i,j)*psi_kanan(i,j);
			}}
			
			element *= dx*dy;
			cout << m << n << " " <<p << q << " " << element << endl;
		}}}}
	
	return 0;
}