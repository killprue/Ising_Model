#include <iostream>
#include <cstdlib>
#include <fstream>
#include <time.h>
using namespace std;

const double gmu0 = 1.26*(9.274009994E-24);
const double Tc = 1040.15;
const double kb = 1.38064852E-23;
const double J = (0.55)*(kb*Tc*log(1+sqrt(2))/2);
const int N = 20;
const int T_range = 150;
const int BigNumber = 500000;

int main() {
    double B = 0.1;
    srand(time(NULL));

    double*** S = new double**[(N+2) * (N+2) * (N+2)];

    for (int i = 0; i < (N+2); i++) {
        S[i] = new double*[N+2];
        for(int j = 0; j < (N+2); j++){
            S[i][j] = new double [N+2];
            for(int k = 0; k < (N+2); k++){
                S[i][j][k] = 0;
            }
        }
    }

    for (int i = 1; i < (N+1); i++) {
        for(int j = 1; j < (N+1); j++){
            for(int k = 1; k < (N+1); k++){
            S[i][j][k] = 1;
            }
        }
    }

    double *Eavg = new double[T_range];
    double *E2avg = new double[T_range];
    double *Mavg = new double[T_range];
    double *M2avg = new double[T_range];
    for(int i = 0; i < T_range; i++){
        Eavg[i] = 0;
        E2avg[i] = 0;
        Mavg[i] = 0;
        M2avg[i] = 0;
    }
    int T = 0;
    for(;T < T_range; T++){
        double E = 0;
        for(int i = 1; i < (N+1); i++) {
            for(int j = 1; j < (N+1); j++){
                for(int k = 1; k < (N+1); k++){
                    E = E - gmu0*S[i][j][k]*B;
                }
            }
        }

        for(int i = 1; i < (N+1); i++) {
            for(int j = 1; j < (N+1); j++){
                for(int k = 1; k < (N+1); k++){
                    E = E - J*S[i][j][k]*(S[i+1][j][k] + S[i-1][j][k] + S[i][j+1][k] + S[i][j-1][k] + S[i][j][k-1] + S[i][j][k+1]);
                }
            }
        }

        double M = 0;
        for(int i = 1; i < (N+1); i++) {
            for(int j = 1; j < (N+1); j++){
                for(int k = 1; k < (N+1); k++){
                    M = M + S[i][j][k];
                }
            }
        }

        double E2 = E*E, M2 = M*M;
        double Esum = E, E2sum = E2, Msum = M, M2sum = M2;
        double kt = kb*T*10;

        for (int iter = 0; iter < BigNumber; iter++){
            int base = 1;
            int range = N;
            int i = base + (rand() % range), j = base + (rand() % range), k = base + (rand() % range);
            S[i][j][k] = -S[i][j][k];
            double DeltaE = -2*gmu0*S[i][j][k]*B - 2*J*S[i][j][k]*(S[i+1][j][k] + S[i-1][j][k] + S[i][j+1][k] + S[i][j-1][k] + S[i][j][k+1] + S[i][j][k-1]);
            double DeltaM = 2 * S[i][j][k];
            if(DeltaE < 0){
                E = E+DeltaE;
                M = M+DeltaM;
                Esum = E+Esum;
                E2sum = E2sum + E*E;
                Msum = Msum + M;
                M2sum = M2sum + M*M;
            }
            else{
                double Prob = ((double) rand() / (RAND_MAX));
                if(Prob < exp(-DeltaE/kt)){
                    E=E+DeltaE;
                    M=M+DeltaM;
                    Esum=Esum+E;
                    E2sum=E2sum+pow(E,2);
                    Msum=Msum+M;
                    M2sum=M2sum+pow(M,2);
                }
                else{
                    S[i][j][k] = -S[i][j][k];
                    Esum = E + Esum;
                    E2sum = E2sum + pow(E,2);
                    Msum=Msum+M;
                    M2sum=M2sum+pow(M,2);

                }
            }
        }
        Eavg[T]=Esum/BigNumber;
        Mavg[T]=Msum/BigNumber;
        E2avg[T]=E2sum/BigNumber;
        M2avg[T]=M2sum/BigNumber;
    }

    ofstream file;
    file.open("isingModel(3d).csv");
    file << "Magnetization,Energy,Magnetization^2,Energy^2" << endl;
    for(int i = 0; i < T_range; i++){
        file << gmu0*Mavg[i] << "," << Eavg[i] << "," << M2avg[i]*gmu0*gmu0 << "," << E2avg[i] << endl;
    }
    file.close();

    for (int i = 0; i < N+2; i++) {
        for (int j = 0; j < N+2; j++) {
            delete[] S[i][j];
        }
    }

    for (int i = 0; i < N+2; i++) {
        delete[] S[i];
    }
    delete[] S;
    delete[] Eavg;
    delete[] E2avg;
    delete[] Mavg;
    delete[] M2avg;

    
}
