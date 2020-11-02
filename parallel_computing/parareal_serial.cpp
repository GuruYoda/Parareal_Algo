#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>

using namespace std;

// Define the problem setup for the differential equation being solved
//Start//

#define A -0.2
#define m 1
#define NUM_THREADS 100
long double T_MIN = 0.0;
long double T_MAX = 10000.0;      
int N_COARSE  = 100;               // (T_max - T_min)*(NUM_THREADS * A)/2
int SCALE_MESH  = 30;
long double Y0  = 10.0;
//End//

// Function to return a 1-D array (Dynamic array  generation method)
double *d_solve(int n_, double t[], double y0_)
{
    double dt = t[1] - t[0];
    double *y = new double[n_]; // solution array

    // Loop to initialize the array
    for (int i = 0; i < n_; i++)
    {
        y[i] = 0; // Initialize the array with zeros
    }

    // Defining the initial value of the system
    y[0] = y0_;

    /*
    // Explicit Euler Difference scheme to solve the IVP
    for (int i = 0; i < n_ - 1; i++)
    {
        y[i + 1] = (1 + dt * A) * y[i];
    }
    */

    // Runge-Kutta Method to solve the IVP
    for (int i = 0 ; i < n_ - 1 ; i++)
    {
        double k1 = A * y[i];
	    double k2 = A * (y[i] + dt * k1/2 );
	    double k = 0.5 * k1 + 0.5 * k2;
	    y[i+1] = y[i] + dt * k;  
    }
    sleep(1);

    return y;
}

void writeArrayTo(ofstream myFile, double arrayToCopy[], int n_len)
{
    if (myFile.is_open())
    {
        for (int i = 0; i < n_len; i++)
        {
            myFile << arrayToCopy[i] << " ";
            myFile << "\n";
        }
    }
}

int main()
{   
    //Start of timer
    double start_time = omp_get_wtime();

    // Time domain over which the differential equation is the be solved
    double t_min = T_MIN;
    double t_max = T_MAX;
    int scale_mesh = SCALE_MESH; // scale size from coarse to fine

    // Initial value
    double y0[1];
    y0[0] = Y0;

    double coarse_start = omp_get_wtime();
    // Corase grid solution to be computed in serial
    int n_sub = NUM_THREADS;
    int n_coarse = N_COARSE;
    //int n_coarse = int(C_DT * (t_max - t_min)/n_sub);
    int n_c = n_sub * n_coarse + 1;
    double dt_c = (t_max - t_min) / (n_c - 1);

    // Discretize time domain with coarse grid
    double t_c[n_c];
    t_c[0] = t_min;
    for (int i = 1; i < n_c; i++)
    {
        t_c[i] = t_min + dt_c * i;
    }

    // Solve IVP on the coarse grid
    int n_ = sizeof(t_c) / sizeof(t_c[0]);
    double *y_c = d_solve(n_, t_c, y0[0]);

    double coarse_end = omp_get_wtime();

    double coarse_time = coarse_end - coarse_start;

    double s_mat[n_sub + 1];    // corrected solution initialized
    double s_mat_prev[n_sub + 1];   // previous coarse grid solution
    double s_mat_new[n_sub + 1];    // current coarse grid solution

    s_mat[0] = y_c[0];
    s_mat_prev[0] = y_c[0];
    s_mat_new[0] = y_c[0];
    for (int i = 1; i < n_sub; i++)
    {
        s_mat[i] = y_c[i * n_coarse];
        s_mat_prev[i] = y_c[i * n_coarse];
        s_mat_new[i] = y_c[i * n_coarse];
    }
    s_mat[n_sub] = y_c[n_];
    s_mat_prev[n_sub] = y_c[n_];
    s_mat_new[n_sub] = y_c[n_];

    //Fine grid solution to be computed in parallel
    int n_fine = scale_mesh * n_coarse;
    int n_f = n_sub * n_fine + 1;
    double dt = (t_c[n_coarse] - t_c[0]) / n_fine;
    //cout << "dt : " << dt << "\n";
    double t[n_f];

    t[0] = t_c[0];
    for (int i = 0; i < n_sub; i++)
    {
        for (int j = (i)*n_fine + 1; j < (i + 1) * n_fine + 1; j++)
        {
            t[j] = t_min + dt * j;
        }
    }

    double yExact[n_f];

    double y[n_f];
    for (int i = 0; i < n_f; i++)
    {
        y[i] = 0;
        yExact[i] = y0[0] * exp(A * t[i]);
    }

    // Defining the parameters for iterations
    double tol = 1e-5;
    double err = 10 * tol;
    int iter = 1;
    int max_iter = 10000;
    double err_vec[max_iter];
    for (int i = 0; i < max_iter; i++)
    {
        err_vec[i] = 0;
    }
    
//    int ind;
//    int start_ind;
    // Iterating over fine grid with better initial guess
    while (err > tol and iter < max_iter)
    {
        int i;         //,t_sub_c,y_sub_c,t_sub,y_sub

        //double fine_start = new double;
        //fine_start = omp_get_wtime();
        for (int i = 0; i < n_sub; i++)
        {
            // fine integrator

            int start_ind = i * n_fine;
            int end_ind = (i + 1) * n_fine;
            double t_sub[n_fine + 1];
            for (int j = 0; j < n_fine + 1; j++)
            {
                t_sub[j] = t[start_ind + j];
            }

            int n_ = sizeof(t_sub) / sizeof(t_sub[0]);
            double *y_sub = d_solve(n_, t_sub, s_mat_prev[i]);

            for (int j = 0; j < n_fine + 1; j++)
            {
                y[start_ind + j] = y_sub[j];
            }
        }
        //double fine_end = new double;
        //fine_end = omp_get_wtime();

        // Update the initial values from fine grid
        for (int i = 1; i < n_sub; i++)
        {
            s_mat[i] = y[i * n_fine];
        }
        s_mat[n_sub] = y[n_f];

        // prediction step
        for (int i = 0; i < n_sub; i++)
        {
            // coarse integrator
            double t_sub_c[n_coarse + 1];
            int start_ind = i*n_coarse;
            for (int j = 0; j < n_coarse + 1; j++)
            {
                t_sub_c[j] = t_c[start_ind + j];
            }

            double *y_sub_c = d_solve(n_coarse+1, t_sub_c, s_mat[i]);
            s_mat_new[i+1] = y_sub_c[n_coarse];
        }

        // Correction step
        for (int i = 1; i < n_sub; i++)
        {
            s_mat[i] = s_mat[i] + s_mat_new[i] - s_mat_prev[i];
            s_mat_prev[i] = s_mat_new[i];
        }
        s_mat[n_sub] = s_mat[n_sub] + s_mat_new[n_sub] - s_mat_prev[n_sub];
        s_mat_prev[n_sub] = s_mat_new[n_sub];

        err = 0.0;
        for (int i = 1; i < n_sub + 1; i++)
        {
            int ind = i * n_fine;
            //if (y[ind] < 10e-5){y[ind] = 0.0;}
	    	//if (s_mat[ind] < 10e-5){s_mat[ind] = 0.0;}
	    	err = err + pow(y[ind] - s_mat[i], 2);
        }

        err = sqrt(dt*err);

        err_vec[iter - 1] = err;
        iter++;
    }
    
    double end_time = omp_get_wtime();
    double full_time = end_time - start_time;
    //double fine_time = fine_end - fine_start;

    string fileName = "solver_summary.txt";
    ofstream myFile;
    myFile.open(fileName, ios::out);
    if (myFile.is_open())
    {
        myFile << "Following is the overall solver summary" << "\n";
        myFile << "A             : " << A << "\n";
        myFile << "m             : " << m << "\n";
        myFile << "T_MIN         : " << T_MIN << "\n";
        myFile << "T_MAX         : " << T_MAX << "\n";
        myFile << "Initial Value : " << Y0 << "\n";
        myFile << "n_coarse      : " << n_coarse << " \n";
        myFile << "Scale Mesh    : " << SCALE_MESH << " \n";
        myFile << "n_fine        : " << n_fine << " \n";
        myFile << "n_sub         : " << n_sub << "\n" ;
        myFile << "dt_coarse     : " << dt_c << "\n";
        myFile << "dt_fine       : " << dt << "\n";
        myFile << "Tolerance     : " << tol << "\n";
        myFile << "Error at exit : " << err << "\n";
        myFile << "Iterations    : " << iter - 1 << "\n";
        myFile << "Time serial   : " << full_time << "\n";
        myFile << "Time coarse   : " << coarse_time << "\n";
            
    }
    myFile.close();

    fileName = "t_coarse.txt";
    myFile.open(fileName, ios::out);
    

    //writeArrayTo(myFile, t_c, n_c);
    if (myFile.is_open())
    {
        for (int i = 0; i < n_c; i++)
        {
            myFile << t_c[i] << "\n";
        }
    }
    myFile.close();

    fileName = "y_coarse.txt";
    myFile.open(fileName, ios::out);
    //writeArrayTo(myFile, y_c, n_c);
    if (myFile.is_open())
    {
        for (int i = 0; i < n_c; i++)
        {
            myFile << y_c[i] << "\n";
        }
    }
    myFile.close();

    fileName = "t_fine.txt";
    myFile.open(fileName, ios::out);
    //writeArrayTo(myFile, t, n_f);
    if (myFile.is_open())
    {
        for (int i = 0; i < n_f; i++)
        {
            myFile << t[i] << "\n";
        }
    }
    myFile.close();

    fileName = "y_fine.txt";
    myFile.open(fileName, ios::out);
    //writeArrayTo(myFile, y, n_f);
    if (myFile.is_open())
    {
        for (int i = 0; i < n_f; i++)
        {
            myFile << y[i] << "\n";
        }
    }
    myFile.close();

    fileName = "y_exact.txt";
    myFile.open(fileName, ios::out);
    //writeArrayTo(myFile, y, n_f);
    if (myFile.is_open())
    {
        for (int i = 0; i < n_f; i++)
        {
            myFile << yExact[i] << "\n";
        }
    }
    myFile.close();

    return 0;
}
