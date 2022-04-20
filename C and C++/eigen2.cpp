#include<cmath>
#include<iostream>
#include<Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

int main(){
	
	// PARAMETER
	int i, j, m, n, p, q, ii, jj;
	int NX = 100, NY = 100, LX = 10, LY = 10;
	float dx = float(LX)/float(NX), dy = float(LY)/float(NY);
	
	float x0, y0;
	x0 = float(LX)/2.0;
	y0 = float(LY)/2.0;
	MatrixXd xo, yo;
	xo = MatrixXd::Constant(NX+1, NY+1, x0);
	yo = MatrixXd::Constant(NX+1, NY+1, y0);
	
	double a = 2/sqrt(LX*LY);
	double pix = M_PI/LX, piy = M_PI/LY;
	
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
	
	//Array tempat menyimpan nilai
	MatrixXd psi_kiri(NX+1, NY+1), psi_kanan(NX+1, NY+1);
	MatrixXd d2psix(NX+1, NY+1), d2psiy(NX+1, NY+1);
	
	// Indexing
	i = 0;
	MatrixXd ind(100, 2);
	for(m = 1; m <= 10; m++){
	for(n = 1; n <= 10; n++){
		ind(i, 0) = m;
		ind(i, 1) = n;
		i++;
	}}	
	
	// tempat menyimpan nilai
	MatrixXd rep_matriks(NX, NY);
	double element;
	
	for(ii=0; ii<NX; ii++){
		m = ind(ii, 0);
		n = ind(ii, 1);
		
		for(i=0; i<NX+1; i++){
		for(j=0; j<NY+1; j++){
			psi_kiri(i, j) = a * sin(m*pix*x(i,j)) * sin(n*piy*y(i,j));
		}}
		
		for(jj=0; jj<NY; jj++){					
			p = ind(jj, 0);
			q = ind(jj, 1);			
			
			for(i=0; i<NX+1; i++){
			for(j=0; j<NY+1; j++){
				psi_kanan(i,j) = a * sin(p*pix*x(i,j)) * sin(q*piy*y(i,j));
			}}
			
			for(i=0; i<NX+1; i++){
			for(j=0; j<NY+1; j++){
				d2psix(i,j) = -0.5 * -(p*piy)*(p*piy) * a * sin(p*piy*x(i,j)) * sin(q*piy*y(i,j));
				d2psiy(i,j) = -0.5 * -(q*piy)*(q*piy) * a * sin(p*piy*x(i,j)) * sin(q*piy*y(i,j));
			}}
			
			element = 0.0;
			for(i=1; i<NX; i++){
			for(j=1; j<NY; j++){
				element += psi_kiri(i,j)*(d2psix(i,j)+d2psiy(i,j)) + psi_kiri(i,j)*V(i,j)*psi_kanan(i,j);
			}}
			
			element *= dx*dy;
			
			//cout << m << n << " " << p << q << " " << element << endl;
			
			rep_matriks(ii, jj) = element;
		}
	}
	
	// cout << ind;
	// cout << rep_matriks;
	EigenSolver<MatrixXd> es(rep_matriks);
	
	cout << "The eigenvalues of Matriks are:" << endl << es.eigenvalues() << endl;
	cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;
	
	return 0;
}