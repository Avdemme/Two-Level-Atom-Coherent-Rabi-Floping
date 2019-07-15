/*This code simulates a two-level atom driven on resonance. It measures the 
 *probability that the atom is in the excited state given that it begins in 
 * the ground start and averages over a given number simulations. We simulate 
 * this atom using the Quantum Monte Carlo Wave Function algorithm as described 
 * in Klaus Molmer's paper "Density Matrices and the Quantum Monte-Carlo Method 
 * in Quantum Optics". The data is a reproduction of figure 3 in this paper. 
 * All matrices were derived analytically.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <complex>

#include "math_fns.h"


using namespace std;

int main(int argc, char** argv) {
    /*First we will initialize our wave-functions:real and imaginary component*/
    std::complex<double> psi_i[2]={{0,0},{1,0}};
    std::complex<double> psi_f[2]={{0,0},{1,0}};
    std::complex<double> psi_n[2]={{0,0},{1,0}};
    
    /*Then we adjust step sizes*/
    double dt = 0.0004;
    double t = 4;
    int step = t/dt;
    int fns = 100;
    
    /*Then we create arrays to store our probabilities
     *and also an array to represent a time scale*/
    double p_rec[fns][step] = {0};
    double t_rec[step] = {0};
    
    /*We create our time axis*/
    for(int i=0; i<step; i++){
        t_rec[i] = i*dt;
    }
    
    /*We initialize our probability*/
    double dp = 0;
    double p = 0;
    
    /*and we set values for our constants*/
    double hbar = 1;
    double gamma = 1;
              
    /*From there we create the real and imaginary parts of our time evolution
     *unitary which we get by taylor expanding our matrix exponential since we
     *were given the Hamiltonian matrix*/
    std::complex<double> U[2][2]={{{1-(gamma*dt),0},{0,-3*gamma*dt}},
                                  {{0,-3*gamma*dt},{1,0}}};
    
    /*We then create our jump operator*/
    std::complex<double> L[2][2] = {{{0,0},{0,0}},{{2*sqrt(gamma),0},{0,0}}};
    
    /*and we must also initialize our random number and its object and counter*/
    double eps = 0;
    double rngind = 0;
    math_fns rngnum;
    
    /*and initialize our norm number and object*/
    double nor = 0;
    math_fns size;        
    
    /*Now we must create multiple wave functions*/
    for(int i=0; i<fns; i++){
        
        /*We reset for every trajectory*/
        psi_i[0] = 0;
        psi_i[1] = 1;
        psi_f[0] = 0;
        psi_f[1] = 1;
        psi_n[0] = 0;
        psi_n[1] = 1;
        dp = 0;
        p = 0;
        eps = 0;
        nor = 0;
        
        /*Now we can implement our QMCWF algorithm to obtain a quantum
         *trajectory*/
        for(int j=0; j<step; j++){
            
            /*We find a random number for this step*/
            eps = rngnum.rng(rngind);
            
            /*Now we must make it so that we start with the wavefunction that
             *was created at the end of the last step*/
            for(int k=0; k<2; k++){
                psi_i[k] = psi_n[k];
            }
            
            /*Now we calculate our dp*/
            dp = gamma*dt*pow(abs(psi_i[0]),2);
            
            /*Now we compare our random number to dp and decide how to evolve
             *our wavefunction*/   
            if(eps <= dp){
                for(int k=0; k<2; k++){
                    psi_f[k] = 0;
                    for(int l=0; l<2; l++){
                        psi_f[k] += L[k][l]*psi_i[l];
                    }
                }
            }        
            else if(eps > dp){
                for(int k=0; k<2; k++){
                    psi_f[k] = 0;
                    for(int l=0; l<2; l++){
                        psi_f[k] += U[k][l]*psi_i[l];
                    }
                }
            }
            
            /*We must normalize our final state*/
            nor = size.norm(psi_f);
            for(int k=0; k<2; k++){
                psi_n[k] = psi_f[k]/nor;
            }
            
            /*We calculate the probability of finding the normalized 
             *state in the excited state*/
            p = pow(abs(psi_n[0]),2);
            
            /*We store the probability in an array to use later*/
            p_rec[i][j] = p;
             
            /*and we advance our random number index*/
            rngind++;    
        }   
    }
    
    /*Now we must find the average probability for each step over all 
     *wavefunctions*/
    double p_avg[step] = {0};
    for(int i=0; i<step; i++){
        for(int j=0; j<fns; j++){
            p_avg[i] += (1./fns)*p_rec[j][i]; 
        }
    }
    
    
    
    
   /*We will create a data file for us to store the data in.*/
   std::remove("QO_data.txt");
   ofstream data;
   data.open ("QO_data.txt", ios::app);
   data.precision(6);
   for(int i=0; i<step; i++){
       data << t_rec[i] << "\t" << p_avg[i] << endl;
   }
   data.close();
    

    

    return 0;
}

