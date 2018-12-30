/* Frank Buono
Calculus of variations program:
Calculates and minimizes area of a shape
prints coordinates to file for visualization
*/

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include "shapetype.h"
#include "grid.h"
using namespace std;

//Declare/Initialize Objects
double RADIUS = 1.0, HEIGHT = 2.0;
ShapeType shape(RADIUS,HEIGHT);
Grid grid(50,30);
int num_unk  = 2 * grid.GetNumberofVertices();
int num_tri  = grid.GetTriangle();
int num_vert = grid.GetNumberofVertices();
int n_grid   = grid.GetN(); //number of horizontal segments
int m_grid   = grid.GetM(); //number of vertical segments
double dz    = shape.GetHeight()/grid.GetN();
int MAX_DEGREE = 14;


//Edge Adjacency Matrix
struct EdgeMatrix
{
        int* aa;
        int* aj;
        int* ai;
};

//calculation constants
double pi   = M_PI;
double dt   = 1.0e-0; //times stepping parameter
int MAX_TSTEP = 200; 
//Dynamically Allocate Arrays
double* coordinatesArr = (double*) malloc(sizeof(double) *
                                         num_unk);

double* gradientArr = (double*) malloc(sizeof(double) *
                                            num_unk);

double* cs_area = (double *) malloc (sizeof(double) * (n_grid + 1));

double * Grad_cs_area = (double *) malloc(sizeof(double)*num_unk);

double * aE  = (double*) malloc( sizeof(double) * 6 * num_tri);
int    * jE  = (int*) malloc( sizeof(int) * 6 * num_tri);
int    * iE  = (int*) malloc( sizeof(int) * (num_tri+1));

double * aET = (double*) malloc( sizeof(double) * 6 * num_tri);
int    * jET = (int*) malloc( sizeof(int) * 6 * num_tri);
int    * iET = (int*) malloc( sizeof(int) * (num_unk+1));

double * aHessPer = (double*) malloc( sizeof(double) * num_unk * MAX_DEGREE);
double * aHessGi  = (double*) malloc( sizeof(double) * num_unk * MAX_DEGREE);
double * aA       = (double*) malloc( sizeof(double) * num_unk * MAX_DEGREE);
int    * jA       = (int*) malloc( sizeof(int) * num_unk * MAX_DEGREE);
int    * iA       = (int*) malloc( sizeof(int) * (num_unk+1));

double * X       = (double*) malloc (sizeof(double) * num_unk);  //coordinate array (the unkown)
double * X0      = (double*) malloc (sizeof(double) * num_unk);  //coordinate array (the unkown)
double   Per;      //perimeter (shape's surface area)
double * GradPer = (double*) malloc (sizeof(double) * num_unk);    //perimeter gradients;
double   Konst   = 1.0e4; //penalty parameter needs dt ~ 1e-3
double   G0      = 313.333;//10.0*10.0*M_PI; //reference cross section area
double * Gi      = (double*) malloc (sizeof(double) * (n_grid+1)); //constraints i = 0,...,n_grid
double * GradGi  = (double*) malloc (sizeof(double) * num_unk);    //lumped constraint gradients
double * V       = (double*) malloc (sizeof(double) * num_unk);    //normalized constraint gradients: V = sqrt(K)GradGi
int    * iV      = (int *) malloc(sizeof(int)*(n_grid+2));         //structuring GradGi into one vector (vs. n_grid + 1 many smaller vectors)

double * F      = (double*) malloc (sizeof(double) * num_unk); //linear system rhs
double * dX     = (double*) malloc (sizeof(double) * num_unk); //linear system solution
double * RR     = (double*) malloc (sizeof(double) * num_unk); //this and following are CG search directions
double * DD      = (double*) malloc (sizeof(double) * num_unk);
double * QQ      = (double*) malloc (sizeof(double) * num_unk);

//Function Headers
void PrintCoordinates(double, double, double, ofstream&);

double triangle_area(double, double, double,
                     double, double, double,
                     double, double, double);

void SetShapeCoordinates(ofstream&);
//Fill Coordinates Array with the coordinate of each vertex
//on the grid
void MinimizeSurfaceArea(ofstream&, int);
//Find the coordinates of surface minimizing shape
//via gradient descent method
//store coordinates in gradient array
//print coordinates to Coordinates file for
//visualization

void TaylorCheck();

void AddEntry(int, int, double,int* ai, int* aj, double* aa);

void Transpose(int m, int n, double *a, int *ia, int *ja, double *b, int *ib, int *jb);

void AtimesX(int m, double * aa, int * aj, int * ai, double * x, double * b);

void ApVVtX(int n, double * aa, int * aj, int * ai, int m, int *iV, double * V, double * x, double * b);

void ATimesB(int m, int n, int p, double *a, int *ia, int *ja, double *b, int *ib, int *jb, double *c, int *ci, int *cj);

double DotProduct(int n, double * x, double *y);

void CGSpecial(int n, int * ia, int *ja, double *aa,
               int m, int *iV, double *V,
               double *x, double *b,
               double *r, double *d, double *q);

void SetCSArea();
//Calculates The crossection at each level N of the Grid
//stores values in cs_area array

void SetTriAreaHessian (double X_in[], double *da, double local_Grad[], double local_Hessian[][9]);

int main(){

//*****************File Setup*****************************
    //open file to print to
    ofstream Vertices_file, Coordinates_file;
    Vertices_file.open("TriangleVert.dat");
    if(Vertices_file.fail())
    {
        cout << "TriangleVert.dat failed to open" << endl;
        exit(0);
    }
    Coordinates_file.open("Coordinates.dat");
    if(Coordinates_file.fail())
    {
        cout << "Coordinates.dat failed to open" << endl;
        exit(0);
    }
    grid.SetTriangleVertices();
//************************************************

    SetShapeCoordinates(Coordinates_file);

    //Assembling edge adjacency matrix

    int i, j, k, jj, t;

    iE[0] = 0;

    for(t = 0; t < num_tri; t++){

        jj = iE[t];

        iE[t+1] = iE[t] + 6;

        i = grid.triangle_vertices[(3 * t)+0];
        j = grid.triangle_vertices[(3 * t)+1];
        k = grid.triangle_vertices[(3 * t)+2];

//      x               y
        jE[jj+0] = 2*i; jE[jj+1] = 2*i + 1;
        jE[jj+2] = 2*j; jE[jj+3] = 2*j + 1;
        jE[jj+4] = 2*k; jE[jj+5] = 2*k + 1;

//      filler; not needed
        aE[jj+0] = 1.0;
        aE[jj+1] = 1.0;
        aE[jj+2] = 1.0;
        aE[jj+3] = 1.0;
        aE[jj+4] = 1.0;
        aE[jj+5] = 1.0;

    }

    Transpose(num_tri, num_unk, aE, iE, jE, aET,iET,jET);

    ATimesB(num_unk, num_tri, num_unk, aET, iET, jET, aE, iE, jE, aA, iA, jA);

    iV[0] = 0;

    for(i = 0; i < n_grid+1; ++i) {

            iV[i+1] = iV[i] + 2*m_grid;

    }

    //done Assembling


    /*

      Start time stepping here

     */
    //    for(i=0;i<n_grid+2;++i) printf("%d\n", iV[i]);
    //    for(i = 0; i < num_unk; i++) printf("%d %4.4g\n",i, X[i]);
    //    exit(1);

    for(int TSTEP = 0; TSTEP <MAX_TSTEP; ++TSTEP) {

      for(i = 0; i < num_unk; i++) X0[i] = X[i];

      MinimizeSurfaceArea(Coordinates_file, TSTEP==0);

    }

    FILE *pfile;
    pfile = fopen("X.dat","w");

    //print results
    for(i=0;i<num_unk/2;++i) {
      double dz = HEIGHT/((double)n_grid);
      int lev = i/m_grid;
      double z = dz*lev;
      fprintf(pfile,"%g %g %g\n", X[2*i], X[2*i+1], z );
    }

    fclose(pfile);

    return 0;
}//end main

void PrintCoordinates(double x_coordinate, double y_coordinate,
                                   double z_coordinate, ofstream& file)
{
        file << x_coordinate << ' ' << y_coordinate
             << ' ' << z_coordinate << endl;
}
double triangle_area(double x0, double y0, double z0,
                     double x1, double y1, double z1,
                     double x2, double y2, double z2 )
//returns the triangular area of the cross product of the
//directional derivatives of the sruface
{

        double dx = (y1-y0)*(z2-z0) - (z1-z0)*(y2-y0);
        double dy = (z1-z0)*(x2-x0) - (x1-x0)*(z2-z0);
        double dz = (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0);

        return (0.5*(sqrt(dx*dx + dy*dy + dz*dz)));

}


void SetShapeCoordinates(ofstream& file)

{
    double da = 0.0;
    double theta;
    double dtheta = 2*pi/grid.GetM();
    double x, y, z, r, dz = HEIGHT/n_grid;
    int row;

    int vertex, i, j;
    double shift = 2.0*RADIUS;

    vertex = 0;
    for(i = 0; i < n_grid+1; ++i) {
      
      z = i*dz;

      for(j = 0; j < m_grid;++j) {

	theta = 2*((double) j) / ((double) (m_grid))*M_PI;

	//shifted circles 
	x = RADIUS*cos(theta) + shift*z/HEIGHT; 
	y = RADIUS*sin(theta);

	//squarish 
        double n_sq = 10;  
	r = RADIUS/pow( pow(cos(theta),n_sq) + pow(sin(theta),n_sq), 1/n_sq);
	
	x = r*cos(theta); y = r*sin(theta); 

	X[2*vertex] = x;  X[2*vertex+1] = y;
	++vertex;
	//	printf("%g %g\n", x, y);
      }      
    }
    /*
    for(int vertex = 0; vertex < num_vert; vertex++)
    {
        theta = grid.GetMCoordinate(vertex) * dtheta;
        row = grid.GetNCoordinate(vertex);

        shape.SetZCoordinate(dz,row);

        x = shape.GetXCoordinate(theta);
        y = shape.GetYCoordinate(theta);
        z = shape.GetZCoordinate();

        X[2*vertex] = x;  X[2*vertex+1] = y;

        PrintCoordinates(x,y,z,file);

    }
    */
}


void MinimizeSurfaceArea(ofstream& Coordinates_file, int START_SIGNAL)
{
    //Placeholders for values stored in the arrays
    double x0,y0,z0; //I ... triangle vertices
    double x1,y1,z1; //J
    double x2,y2,z2; //K
    double xm,ym,dx,dy; //edge midpoints; z-not needed
    int I, J, K, e, i, j, k, l, lev, t, jj, ix, iy, jx, jy, kx, ky, lx, ly;
    int in1[6], in2[6];
    int is_intr;
    int num_tri = grid.GetTriangle();
    double x_coordinate, y_coordinate, z_coordinate;

    //local arrays to do finite differencing.
    double dPer;
    double X_in[9];
    double local_Grad[9];
    double local_Hessian[9][9];

    int MAX_NSTEP = 10;    //max Newton Steps
    double Ntol = 1.0e-5; //Newton tolerance
    double GERR, NERR, NERR0, NTOL = 1.0e-10;
    int NSTEP;
    //use current iterate as initial guess
    for(i = 0; i < num_unk; i++) X[i] = X0[i];


    NSTEP = 0;
    NERR  = NTOL + 1.0;
    
    while(NERR > NTOL) {

    //Assemble values and derivatives relatex to perimeter function
    Per = 0.0; //zero out the perimeter

    for(i = 0; i < num_unk; i++) GradPer[i] = 0.0; //zero out the perimeter gradient

    for(i = 0; i < iA[num_unk]; ++i) aHessPer[i] = 0.0; //zero out the matrix

    for(t = 0; t < num_tri; t++) {

        //get vertices on each triangular area
        I = grid.triangle_vertices[3*t];
        J = grid.triangle_vertices[3*t + 1];
        K = grid.triangle_vertices[3*t + 2];

        //indices of x and y coordinates
        //in coordinates X array
        ix = 2*I; iy = 2*I+1;
        jx = 2*J; jy = 2*J+1;
        kx = 2*K; ky = 2*K+1;

        x0 = X[ix];
        y0 = X[iy];
        z0 = grid.GetNCoordinate(I) * dz;
        x1 = X[jx];
        y1 = X[jy];
        z1 = grid.GetNCoordinate(J) * dz;
        x2 = X[kx];
        y2 = X[ky];
        z2 = grid.GetNCoordinate(K) * dz;

        //pass coordinates of triangle indices to calculate
        //the area of each triangle

        X_in[0] = x0; X_in[1] = y0; X_in[2] = z0;
        X_in[3] = x1; X_in[4] = y1; X_in[5] = z1;
        X_in[6] = x2; X_in[7] = y2; X_in[8] = z2;

        //mapping from local indeces to global indeces
        in1[0] = ix; in1[1] = iy; in1[2] = jx; in1[3] = jy; in1[4] = kx; in1[5] = ky;
        in2[0] = 0;  in2[1] = 1;  in2[2] = 3;  in2[3] = 4;  in2[4] = 6;  in2[5] = 7;

        SetTriAreaHessian (X_in, &dPer, local_Grad, local_Hessian);

        Per += dPer;

        for(i = 0; i < 6; ++i) {

                GradPer[in1[i]] += local_Grad[in2[i]];
                for(j = 0; j < 6; ++j) {

                        AddEntry( in1[i], in1[j], local_Hessian[in2[i]][in2[j]] , iA, jA, aHessPer);
                }

        }

   }

   //Assemble values and derivatives relatex to constraint function(s)

   for(lev = 0; lev < n_grid + 1; lev++) Gi[lev]    = 0.0; //zero out constraints

   for(i = 0; i < num_unk; i++)          GradGi[i]  = 0.0; //zero out the constraint gradient

   for(i = 0; i < iA[num_unk]; ++i)      aHessGi[i] = 0.0; //zero out the constraint hessian

   GERR = 0.0;
   for(lev = 0; lev < n_grid + 1; lev++)
   {
           for(e = 0; e < m_grid; e++)
           {
                   k = (m_grid * lev) + e;
                   l = (k+1)*(e < m_grid-1) + (m_grid*lev)*(e == m_grid-1);
                   kx = 2*k; ky = 2*k+1;
                   lx = 2*l; ly = 2*l+1;

                   //mapping from local indeces to global indeces
                   in1[0] = kx; in1[1] = ky; in1[2] = lx; in1[3] = ly;
                   in2[0] = 0;  in2[1] = 1;  in2[2] = 3;  in2[3] = 4;

                   x0 = X[kx];
                   y0 = X[ky];
                   x1 = X[lx];
                   y1 = X[ly];

                   dx = x1-x0; dy = y1-y0;
                   xm = 0.5*(x0 + x1);
                   ym = 0.5*(y0 + y1);

                   Gi[lev] += 0.5*((xm*dy) - (ym*dx));
                   GradGi[kx] +=  0.5 * y1;
                   GradGi[ky] += -0.5 * x1;
                   GradGi[lx] += -0.5 * y0;
                   GradGi[ly] +=  0.5 * x0;

                   AddEntry( kx, ly,  0.5 , iA, jA, aHessGi);
                   AddEntry( ky, lx, -0.5 , iA, jA, aHessGi);
                   AddEntry( lx, ky, -0.5 , iA, jA, aHessGi);
                   AddEntry( ly, kx,  0.5 , iA, jA, aHessGi);

           }
	   if(START_SIGNAL) G0 = Gi[0];
	   //   printf("%d %g %g\n", lev, Gi[lev], G0);
	   GERR += (Gi[lev]-G0)*(Gi[lev]-G0);
   }
   //exit(1);
   //Finished Assembling Derivative Matrices

   //Right Hand Side/Filling out matrix A/Zeroing out boundary values

   int i_is_intr, j_is_intr;
   for(i = 0; i < num_unk; ++i)
   {
     i_is_intr = (2*m_grid <= i)&&(i < 2*m_grid*n_grid);
     V[i] = sqrt( dt * Konst) * GradGi[i] * i_is_intr;

     lev = i/(2*m_grid);

     F[i] = X0[i] - X[i] - dt*(GradPer[i] + Konst*(Gi[lev] - G0)*GradGi[i]);
     F[i] *= i_is_intr;
     
     for(jj = iA[i]; jj < iA[i+1]; jj++)
       {
	 j = jA[jj];

	 j_is_intr = (2*m_grid <= j)&&(j < 2*m_grid*n_grid);

	 aA[jj] = (i == j) + dt * (aHessPer[jj] + Konst * (Gi[lev] - G0) * aHessGi[jj])*i_is_intr*j_is_intr;
       }

     //     printf("%d %d %2.2g\n", i, i_is_intr, X[i]);

   }

   //Something's not right; running a Taylor expansion check....

   //   TaylorCheck();   exit(1);


   //Solve (A + VVt)dX = F where VVt is actually the tensor sum v1v1t + ... + vmvmt

   CGSpecial(num_unk, iA, jA, aA,
             n_grid+1, iV, V,
             dX, F,
             RR, DD, QQ);

   // validity check
   //ApVVtX(num_unk, aA, jA, iA, n_grid+1, iV, V, dX, RR);
   // for(i=0;i<num_unk;++i) printf("%g\n", fabs(F[i]- RR[i]));
   // exit(1);
   //   for(i=0;i<num_unk;++i) printf("%g\n", V[i]);
   /*
     for(i=0;i<num_unk;++i) {
     for(j=iA[i];j<iA[i+1];++j) {
       printf("%d %d %g\n", i, jA[j], aA[j]);
     }
   }
   */
   //for(i=0;i<num_unk;++i) printf("%g %g\n", dX[i], RR[i]);
   //exit(1);
   
  

   //Update the current solution   
   NERR = 0.0;
   for(i=0;i<num_unk;++i) {
     X[i] += dX[i];
     NERR += dX[i]*dX[i];
     //     printf("%g\n", dX[i]);
   }

   NERR = sqrt(NERR);
   if(NSTEP==0) NERR0 = NERR;

   /*
   if(NSTEP > 7) {
     dt *= 0.5;
     NSTEP = 0;
     NERR = NTOL+1;
     for(i=0;i<num_unk;++i) X[i] = X0[i];
   }

   if(NERR < NTOL & NSTEP < 5) {
     dt *= 2;
   }
   */

   ++NSTEP;

    printf("%d %g %g %g %g %g\n", NSTEP, NERR, NERR0, Per, 0.5*Konst*GERR, Per + 0.5*Konst*GERR);
    }



}


void TaylorCheck() {


  double eps = 1.0e-5;
  
  //Placeholders for values stored in the arrays
  double x0,y0,z0; //I ... triangle vertices
  double x1,y1,z1; //J
  double x2,y2,z2; //K
  double xm,ym,dx,dy; //edge midpoints; z-not needed
  int I, J, K, e, i, j, k, l, lev, t, jj, ix, iy, jx, jy, kx, ky, lx, ly;
  int in1[6], in2[6];
  int is_intr;
  int num_tri = grid.GetTriangle();
  double x_coordinate, y_coordinate, z_coordinate;
  
  //local arrays to do finite differencing.
  double dPer;
  double X_in[9];
  double local_Grad[9];
  double local_Hessian[9][9];

  double Per0 = Per;
  double * Gi0 = (double *) malloc(sizeof(double)*(n_grid+1));
  for(i=0;i<n_grid+1;++i) Gi0[i] = Gi[i];

  srand( 1 );

  for(i=0;i<num_unk;++i) {

    X0[i] = X[i];
    dX[i] = eps*(2.0*((double) rand())/((double) RAND_MAX)-1.0);
    X[i]  = X0[i] + dX[i];

    //    printf("%g\n", dX[i]);
    
  }
  Per = 0.0;
  for(t = 0; t < num_tri; t++) {

    //get vertices on each triangular area
    I = grid.triangle_vertices[3*t];
    J = grid.triangle_vertices[3*t + 1];
    K = grid.triangle_vertices[3*t + 2];
    
    //indices of x and y coordinates
    //in coordinates X array
    ix = 2*I; iy = 2*I+1;
    jx = 2*J; jy = 2*J+1;
    kx = 2*K; ky = 2*K+1;

    x0 = X[ix];
    y0 = X[iy];
    z0 = grid.GetNCoordinate(I) * dz;
    x1 = X[jx];
    y1 = X[jy];
    z1 = grid.GetNCoordinate(J) * dz;
    x2 = X[kx];
    y2 = X[ky];
    z2 = grid.GetNCoordinate(K) * dz;
    
    //pass coordinates of triangle indices to calculate
    //the area of each triangle
    
    X_in[0] = x0; X_in[1] = y0; X_in[2] = z0;
    X_in[3] = x1; X_in[4] = y1; X_in[5] = z1;
    X_in[6] = x2; X_in[7] = y2; X_in[8] = z2;

    //mapping from local indeces to global indeces
    in1[0] = ix; in1[1] = iy; in1[2] = jx; in1[3] = jy; in1[4] = kx; in1[5] = ky;
    in2[0] = 0;  in2[1] = 1;  in2[2] = 3;  in2[3] = 4;  in2[4] = 6;  in2[5] = 7;

    SetTriAreaHessian (X_in, &dPer, local_Grad, local_Hessian);    
    Per += dPer;
    
   }
  
   //Assemble values and derivatives relatex to constraint function(s)

   for(lev = 0; lev < n_grid + 1; lev++) Gi[lev]    = 0.0; //zero out constraints
   
   for(lev = 0; lev < n_grid + 1; lev++)
   {
           for(e = 0; e < m_grid; e++)
           {
                   k = (m_grid * lev) + e;
                   l = (k+1)*(e < m_grid-1) + (m_grid*lev)*(e == m_grid-1);
                   kx = 2*k; ky = 2*k+1;
                   lx = 2*l; ly = 2*l+1;

                   //mapping from local indeces to global indeces
                   in1[0] = kx; in1[1] = ky; in1[2] = lx; in1[3] = ly;
                   in2[0] = 0;  in2[1] = 1;  in2[2] = 3;  in2[3] = 4;

                   x0 = X[kx];
                   y0 = X[ky];
                   x1 = X[lx];
                   y1 = X[ly];

                   dx = x1-x0; dy = y1-y0;
                   xm = 0.5*(x0 + x1);
                   ym = 0.5*(y0 + y1);

                   Gi[lev] += 0.5*((xm*dy) - (ym*dx));

           }

	   //printf("%16.16g, %16.16g\n", Gi[lev],Gi0[lev]);
   }


   double E0 = Per0;
   double E  = Per;
   for(lev=0;lev<n_grid+1;++lev) {

     E0 += 0.5*Konst*(Gi0[lev]-G0)*(Gi0[lev]-G0);
     E  += 0.5*Konst*(Gi[lev] -G0)*(Gi[lev] -G0);

   }

   for(i=0;i<num_unk;++i) {

     lev = i/(2*m_grid);

     E0 += (GradPer[i] + Konst*(Gi0[lev]-G0)*GradGi[i])*dX[i];

   }

   printf("%g %g %g\n", E0, E, (E-E0)/(pow(eps,2)));
   //two terms was enough; we figured it out already...

}

void SetTriAreaHessian (double X_in[], double *dPer, double local_Grad[], double local_Hessian[][9])
{

    double Xi, Xj;
    double h = 1.0e-4;

    int n_var = 9;

    double X[n_var];
    double temp[2][2];
    double diag[3];

    for(int i = 0; i < n_var; i++){
        for(int j = i+1; j < n_var; j++)
        {

            for(int k=0; k<n_var; ++k) X[k] = X_in[k];
            X[i] -= h; X[j] += h;
            temp[0][0] = triangle_area(X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8]);
            X[i] = X_in[i]; X[j] = X_in[j];

            X[i] += h; X[j] += h;
            temp[0][1] = triangle_area(X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8]);
            X[i] = X_in[i]; X[j] = X_in[j];

            X[i] -= h; X[j] -= h;
            temp[1][0] = triangle_area(X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8]);
            X[i] = X_in[i]; X[j] = X_in[j];

            X[i] += h; X[j] -= h;
            temp[1][1] = triangle_area(X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8]);
            X[i] = X_in[i]; X[j] = X_in[j];

            local_Hessian[i][j] = (temp[0][1] - temp[0][0] - temp[1][1] + temp[1][0])/(4.0*h*h);

            local_Hessian[j][i] = local_Hessian[i][j];

        }

        for(int k=0; k<n_var; ++k) X[k] = X_in[k];
        X[i] -= h;
        diag[0] = triangle_area(X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8]);
        X[i] = X_in[i];
        X[i] += 0.0;
        diag[1] = triangle_area(X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8]);
        X[i] = X_in[i];

        X[i] += h;
        diag[2] = triangle_area(X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8]);

        local_Hessian[i][i] = (diag[2] + diag[0] - 2.0*diag[1])/(h*h);

        local_Grad[i] = (diag[2] - diag[0])/(2.0*h);

    }

    *dPer = diag[1];

}

void AddEntry(int I , int J , double val, int* ai, int* aj, double* aa)
{
    for(int j = ai[I]; j < ai[I+1]; j++)
    {
        if( aj[j] == J) aa[j] += val;
    }
}

void Transpose(int m, int n, double *a, int *ia, int *ja, double *b, int *ib, int *jb) {

  int i, j, jj, k, nnz;
  double temp;
  nnz = ia[m];

  //row_ptrs set to zero
  for(i=0;i<n+1;++i) ib[i] = 0;

  //partition the column index, then make row_ptr point to
  for(i=0;i<m;++i) {
    for(j=ia[i];j<ia[i+1];++j) {
      jj = ja[j];
      ib[jj+1] = ib[jj+1] + 1;
    }
  }
  ib[0] = 0;

  // add row_ptr
  for(i=0;i<n;++i) {
    ib[i+1] = ib[i+1] + ib[i];
  }

  // set col_ind to -1
  for(i=0;i<nnz;++i) jb[i]=-1;

  for(i=0;i<m;++i) {
    for(j=ia[i];j<ia[i+1];++j) {
      jj = ja[j];		//partition the col_ind according to original row_ptr?
      temp = a[j];
      for(k=ib[jj];k<ib[jj+1];++k){
	/*printf("%d,%d,%d,%d, %d,%d\n",i,jj,(*ib)[jj],k,(*ib)[jj+1],(*jb)[k]);*/
	if(jb[k]==-1) {
	  jb[k] = i;
	  b[k] = temp;
	  break;
	}
	/*printf("%d,%d,%d,%d,%d,%d\n",i,jj,(*ib)[jj],k,(*ib)[jj+1],(*jb)[k]);*/
      }
    }
  }

}

void AtimesX(int m, double * aa, int * aj, int * ai, double * x, double * b)
{
    int column, j;
    for( int i = 0; i < m; i ++)
    {
        b[i] = 0;
        for(int j  = ai[i]; j < ai[i+1]; j ++)
        {
            column = aj[j];
            b[i] += aa[j] * x[column];
        }
    }
}

void ApVVtX(int n, double * aa, int * aj, int * ai, int m, int *iV, double * V, double * x, double * b)
{

        //A is n by n mantrix, x is vector length n

        //here VVt means v1v1t + v1v1t + ... vmvmt where we store v1,...,vm as
        // V = [v1 | v2 | ... | vm] with pointer iV

        int i, jj, j;
        double dot;


        for(i = 0; i < n; i ++)
        {
            b[i] = 0;
            for(jj  = ai[i]; jj < ai[i+1]; jj ++)
            {
                j = aj[jj];
                b[i] += aa[jj] * x[j];
            }
        }
	
        for(i = 0; i < m; ++i)
        {
                dot = 0.0;
                for(j = iV[i]; j < iV[i+1]; ++j)
                {
                        dot += V[j]*x[j];
			//printf("%d %d %d %g %g %g\n", m, i, j, dot, V[j], b[j]);
		}

                for(j = iV[i]; j < iV[i+1]; ++j)
                {
		  b[j] += V[j]*dot;
                }

        }


}

void ATimesB(int m, int n, int p, double *a, int *ia, int *ja, double *b, int *ib, int *jb, double *c, int *ci, int *cj) {

  //added row size n, so we can also use only upper left parts
  //of matrices

  int maxdim;
  maxdim = (m>=n)*m + (m<n)*n;

  int i, j, jj, k, kk, nnz, l, ll, event;
  for(i = 0; i <= maxdim; ++i) cj[i] = -1;

  /*count total occurences of connections
    ci keeps track of number of occurences per row*/
  nnz = 0;
  for(i = 0; i < m; ++i) {
    for(j = ia[i]; j < ia[i+1]; ++j) {
      jj = ja[j];
      if(jj<n) {
	for(k = ib[jj]; k < ib[jj+1]; ++k){
	  kk = jb[k];
	  if(kk < p) {
	    //kk already seen? then ci[kk] = i
	    if(cj[kk] != i) {
	      ++nnz;
	      cj[kk] = i; //kk seen this turn, don't count again
	    }
	  }//end if
	}
      }//end if
    }
  }

  for(i = 0; i < nnz; ++i) c[i] = 0.0;

  ci[0] = 0;
  //go through again and retrieve products
  for(i = 0; i < m; ++i) {
    ci[i+1] = ci[i];
    for(j = ia[i]; j < ia[i+1]; ++j) {
      jj = ja[j];

      if(jj<n) {
	for(k = ib[jj]; k < ib[jj+1]; ++k){
	  kk = jb[k];

	  if(kk<p) {
	    //kk already seen? then add product
	    event = 0;
	    for(l = ci[i]; l < ci[i+1]; ++l) {
	      ll = cj[l];
	      if(kk == ll) {
		event = 1;
		c[l] = c[l] + a[j]*b[k];
	      }
	    }
	    //kk not seen? then iterate ia and add to ja and add product
	    if(event == 0) {
	      cj[ci[i+1]] = kk;
	      ci[i+1] = ci[i+1]+1;
	      c[l] = c[l] + a[j]*b[k];
	    }
	  } //end if
	}
      }//end if
    }

  }


}

double DotProduct(int n, double * x, double *y)
{
    double product = 0;
    for(int i = 0; i < n ; i++)
    {
        product+= x[i]*y[i];
    }
    return product;
}

void CGSpecial(int n, int * ia, int *ja, double *aa,
               int m, int *iV, double *V,
               double *x, double *b,
               double *r, double *d, double *q)
{

        // Conjugate Gradient data members
        double TOL = 1.0e-25;
        double dnew, dnaught, dold, alpha,beta;
        int STEP, i, j, k, l;
        int STEPMAX = 300;
	ApVVtX(n, aa, ja, ia, n_grid+1, iV, V, x, r);

	  //        AtimesX(n, aa, ja, ia, x, r); //r = Ax

        //r = b - Ax
        //d = r
        for(int i = 0; i < n; i++)
        {
                r[i] = b[i] - r[i];
                d[i] = r[i];
        }


        dnew = DotProduct(n,r,r);
        dnaught = dnew;

        //cout << "dnew before loop: " << dnew << endl;
	STEP = 0;
        while(STEP < STEPMAX && dnew > (TOL*TOL*dnaught))
        {

	  ApVVtX(n, aa, ja, ia, n_grid+1, iV, V, d, q);
	  //AtimesX(n,aa,ja,ia,d,q);

	  alpha = dnew / DotProduct(n,d,q);

	  // x += alpha*d
	  for(int i = 0; i < n; i++)
	    {
	      x[i] += alpha * d[i];
	    }
	  
	  if(!(STEP%50) && STEP != 0){
	    //cout << "STEP: " << STEP << " is divisible by 50\n";
	    ApVVtX(n, aa, ja, ia, n_grid+1, iV, V, x, r);
	    //                        AtimesX(n,aa,ja,ia,x,r); //r = Ax
	    for(int i = 0; i < n; i++)
	      {
		r[i] = b[i] - r[i];
	      }
	  }
	  else{
	    for(int i = 0; i < n; i++)
	      {
		
		r[i] -= alpha * q[i];
	      }
	  }
	  dold = dnew;
	  dnew = DotProduct(n,r,r);
	  beta = dnew/dold;
	  
	  //	  cout << "dnew iteration " << STEP << ": " << dnew << endl;
	  //cout << "epsilonsquareddnaught: " << TOL*TOL*dnaught << endl;
	  
	  for(int i = 0; i < n; i++)
	    {
	      d[i] = r[i] + beta*d[i];
	    }
	  STEP++;
        }


}

void SetCSArea()
{
        //used to acces coordinates array
        int k,l,kx,ky,lx,ly;
        double x0,y0,x1,y1;
        double dx,dy;
        double xm,ym;

        for(int n = 0; n <= n_grid; n++)
        {
                for(int i = 0; i < m_grid; i++)
                {
                        k = (m_grid * n) + i;
                        l = (k+1)*(i < m_grid-1) + (m_grid*n)*(i == m_grid-1);
                        kx = 2*k; ky = 2*k+1;
                     	lx = 2*l; ly = 2*l+1;

                     	x0 = coordinatesArr[kx];
                     	y0 = coordinatesArr[ky];
                     	x1 = coordinatesArr[lx];
                     	y1 = coordinatesArr[ly];

                     	dx = x1-x0; dy = y1-y0;
                     	xm = 0.5*(x0 + x1);
                     	ym = 0.5*(y0 + y1);

                     	cs_area[n] += (xm*dy) - (ym*dx);
                        Grad_cs_area[kx] +=  0.5 * y1;
                        Grad_cs_area[ky] += -0.5 * x1;
                        Grad_cs_area[lx] += -0.5 * y0;
                        Grad_cs_area[ly] +=  0.5 * x0;

                }
                cs_area[n] = cs_area[n]/2.0;
        }
}
