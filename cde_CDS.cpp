#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <setjmp.h>
#include <iomanip>

using namespace std;

// Domain length
double domainLength(vector<double> domainCorners) {
  return domainCorners[1]-domainCorners[0];
}

// Grid size
double gridSize(double length, int N){
  return length/double(N);
}

// Position of control volumes
vector<double> cvPosition(int N, double d, vector<double> domainCorners){
  vector<double> posCV(N+2, 0.0);
  posCV[0] = domainCorners[0];
  posCV[N+1] = domainCorners[1];
  for(int i=1; i<=N; i++){
    posCV[i]=domainCorners[0]+i*d;
  }
  return posCV;
}

// Position of nodes
vector<double> nodePosition(vector<double> posCV, int N, vector<double> domainCorners){
  vector<double> posP(N+2, 0.0);
  posP[0] = domainCorners[0];
  posP[N+1] = domainCorners[1];
   for(int i=1; i<=N+1; i++){
     posP[i]=(posCV[i]+posCV[i-1])/2.0;
   }
  return posP;
}

double velocityU(double X, double Y) {
    double u=2*Y*(1-X*X);
    return u;
}

double velocityV(double X, double Y) {
    double v=-2*X*(1-Y*Y);
    return v;
}

//Numerical Solution of the Convection-Diffusion Equation
int main(){
  double PHI_o=0.0;                 // Initial conditions
  double alpha=10.0;                // Alpha coefficient
  double gamma=10.0;                // Ratio between rho and GAMMA

  int N=100, M=50;                  // Number of control volumes; mesh

  //Time definition
  double dt=1.0e-5;
    double time=0.0;

  vector<double> domainX = {-1, 1}, domainY = {0, 1};

  double lengthX = domainLength(domainX), lengthY = domainLength(domainY);

  double dx=gridSize(lengthX,N), dy=gridSize(lengthY,M);

  vector<double> xCV=cvPosition(N, dx, domainX), yCV=cvPosition(M, dy, domainY);
  vector<double> xP=nodePosition(xCV, N, domainX), yP=nodePosition(yCV, M, domainY);

  //Solver definition
  double timeTolerance=1e-6;        // Steady-state convergence
  double timeChange=0.0;            // Maximum difference between previous and current time step

  // Preallocate coefficients
  double aE[N+2][M+2], aW[N+2][M+2], aS[N+2][M+2], aN[N+2][M+2], bP[N+2][M+2], aP[N+2][M+2];

  // Initialize property Φ vectors
  double PHIo[N+2][M+2];     // Initial property map
    for (int i = 0; i <= N+1; i++) {
        for (int j = 0; j <= M+1; j++) {
            PHIo[i][j] = PHI_o;
        }
    }

    // Boundary conditions for PHI - Stays constant in time
    // Upper wall
    for (int i=1; i<=N;i++){
        PHIo[i][M+1]=1-tanh(alpha);
    }

    // Left wall
    for (int j=1;j<=M+1;j++){
        PHIo[0][j]=1-tanh(alpha);
    }

    // Right wall
    for (int j=1; j<=M+1;j++){
        PHIo[N+1][j]=1-tanh(alpha);
    }

    // Inlet and outlet - Initial conditions; Inlet will remain constant in time
    for (int i=0; i<=N+1;i++) {
        double x=xP[i];
        if (x<=0) {
            PHIo[i][0]=1+tanh((2.0*x+1)*alpha);
        } else {
            PHIo[i][0]=PHIo[i][1];
        }
    }

    double PHIpn[N+2][M+2];       // Initialize property vector at the previous time step t(n) at a node P
    for (int i = 0; i <= N+1; i++) {
        for (int j = 0; j <= M+1; j++) {
            PHIpn[i][j] = PHIo[i][j];
        }
    }

    double PHIpn_1[N+2][M+2];
    for (int i = 0; i <= N+1; i++) {
        for (int j = 0; j <= M+1; j++) {
            PHIpn_1[i][j] = PHIpn[i][j];
        }
    }

    // Inner coefficients
          for (int i=2; i<N;i++){
            for (int j=2; j<M;j++){
                double ue=velocityU(xCV[i],yP[j]);
                double uw=velocityU(xCV[i-1], yP[j]);
                double vn=velocityV(xP[i],yCV[j]);
                double vs=velocityV(xP[i], yCV[j-1]);

                aE[i][j]=-gamma*ue/(2.0*dx)+1.0/(dx*dx);
                aW[i][j]= gamma*uw/(2.0*dx)+1.0/(dx*dx);
                aN[i][j]=-gamma*vn/(2.0*dy)+1.0/(dy*dy);
                aS[i][j]= gamma*vs/(2.0*dy)+1.0/(dy*dy);
                bP[i][j]=(-gamma*ue/(2.0*dx)-2.0/(dx*dx)+gamma*uw/(2.0*dx) -
                         gamma*vn/(2.0*dy)-2.0/(dy*dy)+gamma*vs/(2.0*dy)+gamma/dt);
                aP[i][j]=gamma/dt;
            }
        }

        // Upper boundary nodes
        for (int i=2; i<=N-1;i++){
            double ue=velocityU(xCV[i],yP[M]);
            double uw=velocityU(xCV[i-1], yP[M]);
            double vn=velocityV(xP[i],yCV[M]);
            double vs=velocityV(xP[i], yCV[M-1]);

            aE[i][M]=-gamma*ue/(2.0*dx)+1.0/(dx*dx);
            aW[i][M]= gamma*uw/(2.0*dx)+1.0/(dx*dx);
            aN[i][M]=-gamma*vn/dy+2.0/(dy*dy);
            aS[i][M]= gamma*vs/(2.0*dy)+1.0/(dy*dy);
            bP[i][M]=(-gamma*ue/(2.0*dx)-2.0/(dx*dx)+gamma*uw/(2.0*dx) -
                     3.0/(dy*dy)+gamma*vs/(2.0*dy)+gamma/dt);
            aP[i][M]=gamma/dt;
          }

        // Left boundary nodes
          for (int j=2; j<=M;j++) {
              double ue=velocityU(xCV[1],yP[j]);
              double uw=velocityU(xCV[0], yP[j]);
              double vn=velocityV(xP[1],yCV[j]);
              double vs=velocityV(xP[1], yCV[j-1]);

              aE[1][j]=-gamma*ue/(2.0*dx)+1.0/(dx*dx);
              aS[1][j]= gamma*vs/(2.0*dy)+1.0/(dy*dy);
              aW[1][j]= gamma*uw/dx+2.0/(dx*dx);
              aP[1][j]=gamma/dt;
              if (j==M) {
                  aN[1][j]=-gamma*vn/dy+2.0/(dy*dy);
                  bP[1][j]=(-gamma*ue/(2.0*dx)-3.0/(dy*dy) -
                       3.0/(dx*dx)+gamma*vs/(2.0*dy)+gamma/dt);
              } else {
                  aN[1][j]=-gamma*vn/(2.0*dy)+1.0/(dy*dy);
                  bP[1][j]=(-gamma*ue/(2.0*dx)-2.0/(dy*dy)-gamma*vn/(2.0*dy) -
                       3.0/(dx*dx)+gamma*vs/(2.0*dy)+gamma/dt);
              }
            }

      // Right boundary nodes
      for (int j=2; j<=M;j++) {
          double ue=velocityU(xCV[N],yP[j]);
          double uw=velocityU(xCV[N-1], yP[j]);
          double vn=velocityV(xP[N],yCV[j]);
          double vs=velocityV(xP[N], yCV[j-1]);

          aE[N][j]=-gamma*ue/dx+2.0/(dx*dx);
          aS[N][j]= gamma*vs/(2.0*dy)+1.0/(dy*dy);
          aW[N][j]= gamma*uw/(2.0*dx)+1.0/(dx*dx);
          aP[N][j]=gamma/dt;
          if (j==M) {
              aN[N][j]=-gamma*vn/dy+2.0/(dy*dy);
              bP[N][j]=(gamma*uw/(2.0*dx)-3.0/(dy*dy) -
                   3.0/(dx*dx)+gamma*vs/(2.0*dy)+gamma/dt);
          } else {
              aN[N][j]=-gamma*vn/(2.0*dy)+1.0/(dy*dy);
              bP[N][j]=(gamma*uw/(2.0*dx)-2.0/(dy*dy)-gamma*vn/(2.0*dy) -
                   3.0/(dx*dx)+gamma*vs/(2.0*dy)+gamma/dt);
          }
      }

        // Bottom boundary nodes
        for (int i=1;i<=N;i++) {

            double ue=velocityU(xCV[i],yP[1]);
            double uw=velocityU(xCV[i-1], yP[1]);
            double vn=velocityV(xP[i],yCV[1]);
            double vs=velocityV(xP[i], yCV[0]);

            double x=xP[i];
            if (x<=0) {
                aE[i][1]=-gamma*ue/(2.0*dx)+1.0/(dx*dx);
                aN[i][1]=-gamma*vn/(2.0*dy)+1.0/(dy*dy);
                aS[i][1]= gamma*vs/dy+2.0/(dy*dy);
                aP[i][1]=gamma/dt;
                if (i==1){
                    aW[i][1]= gamma*uw/dx+2.0/(dx*dx);
                    bP[i][1]=(-gamma*ue/(2.0*dx)-3.0/(dx*dx) -
                         3.0/(dy*dy)-gamma*vn/(2.0*dy)+gamma/dt);
                } else {
                    aW[i][1]= gamma*uw/(2.0*dx)+1.0/(dx*dx);
                    bP[i][1]=(-gamma*ue/(2.0*dx)-2.0/(dx*dx)+gamma*uw/(2.0*dx) -
                         3.0/(dy*dy)-gamma*vn/(2.0*dy)+gamma/dt);
            }
            } else {
                aW[i][1]= gamma*uw/(2.0*dx)+1.0/(dx*dx);
                aN[i][1]=-gamma*vn/(2.0*dy)+1.0/(dy*dy);
                aS[i][1]= 0.0;
                aP[i][1]=gamma/dt;
                if (i==N) {
                    aE[i][1]=-gamma*ue/dx+2.0/(dx*dx);
                    bP[i][1]=(gamma*vs/dy-3.0/(dx*dx)+gamma*uw/(2.0*dx) -
                         1.0/(dy*dy)-gamma*vn/(2.0*dy)+gamma/dt);
                } else {
                    aE[i][1]=-gamma*ue/(2.0*dx)+1.0/(dx*dx);
                    bP[i][1]=(-gamma*ue/(2.0*dx)-2.0/(dx*dx)+gamma*uw/(2.0*dx) -
                         1.0/(dy*dy)-gamma*vn/(2.0*dy)+gamma*vs/dy+gamma/dt);
                }
            }
        }

  // Start the time loop
  do {
      time+=dt;
      cout<<"TIME = "<<time<<endl;
        timeChange=0.0;

            // Internal nodes
            for (int i=1;i<=N;i++){
                for (int j=1;j<=M;j++){
                    PHIpn_1[i][j]=(aE[i][j]*PHIpn[i+1][j]+aW[i][j]*PHIpn[i-1][j]+
                        aN[i][j]*PHIpn[i][j+1]+aS[i][j]*PHIpn[i][j-1]+PHIpn[i][j]*bP[i][j])/aP[i][j];
                }
            }

      // Set new values for inlet/outlet
      // Inlet and outlet - Initial conditions; Inlet will remain constant in time
      for (int i=0; i<=N+1;i++) {
          double x=xP[i];
          if (x<=0) {
              continue;
          } else {
              PHIpn_1[i][0]=PHIpn_1[i][1];
          }
      }

        for (int i = 1; i < N ; i++) {
          for (int j = 1; j < M;  j++) {
              timeChange = max(timeChange, fabs(PHIpn_1[i][j]-PHIpn[i][j]));
          }
        }

      for (int i = 0; i <= N+1; i++) {
          for (int j = 0; j <= M+1; j++) {
              PHIpn[i][j] = PHIpn_1[i][j];
          }
      }
  } while(timeChange>=timeTolerance);

    // Plotting data
    ofstream phiMap("PHI_data.csv");

    if (phiMap.is_open()) {
        for (int j = 0; j <= M + 1; j++) {
            for (int i = 0; i <= N + 1; i++) {
                phiMap << PHIpn[i][j];
                if (i < N + 1) phiMap << ",";
            }
            phiMap << "\n";
        }
        phiMap.close();
        cout << "Property data saved to 'PHI_data.csv'\n";
    } else {
        cerr << "Error: Unable to open file for writing\n";
    }

    ofstream phiOutlet("PHI_outlet.csv");
    if (phiOutlet.is_open()) {
        for (int i = 0; i <= N + 1; i++) {
            double x=xP[i];
            phiOutlet << x << "," << PHIpn[i][0] << "\n";
        }
        phiOutlet.close();
    }

    vector<double> target_x = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    ofstream phiInterpolated("PHI_interpolated.txt");
    phiInterpolated << std::fixed << std::setprecision(3); // Set decimal precision

    if (phiInterpolated.is_open()) {
        for (double x_target : target_x) {
            int i_left = 0;
            while (i_left < N && xP[i_left] < x_target) {
                i_left++;
            }

            if (i_left == 0) {
                phiInterpolated << x_target << " " << PHIpn[i_left][0] << "\n";
                continue;
            }

            int i_right = i_left;
            i_left--; // Now xP[i_left] < x_target < xP[i_right]

            double x1 = xP[i_left], x2 = xP[i_right];
            double phi1 = PHIpn[i_left][0], phi2 = PHIpn[i_right][0];

            double phi_interp = phi1 + (phi2 - phi1) * (x_target - x1) / (x2 - x1);

            phiInterpolated << x_target << " " << phi_interp << "\n";
        }
        phiInterpolated.close();
        cout << "Interpolated PHI values saved to 'PHI_interpolated.txt'\n";
    } else {
        cerr << "Error: Unable to open file for writing\n";
    }

  return 0;
}


