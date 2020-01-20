#include <iostream>
#include <fstream>
#include <string>
#include "sac.hpp"

using namespace std;

void Mat_size(MatrixXd m){
    cout << "rows = " << m.rows() << " cols == "<< m.cols() << endl;
}

void Vec_size(VectorXd v){
    cout << "size = " << v.size() << endl;
}

MatrixXd readMat(string file, int rows, int cols) {

  ifstream in(file);
  
  string line;

  int row = 0;
  int col = 0;

  MatrixXd res = MatrixXd(rows, cols);

  if (in.is_open()) {

    while (getline(in, line)) {

      char *ptr = (char *) line.c_str();
      int len = line.length();

      col = 0;

      char *start = ptr;
      for (int i = 0; i < len; i++) {

        if (ptr[i] == ' ') {
          res(row, col++) = atof(start);
          start = ptr + i + 1;
        }
      }
      res(row, col) = atof(start);

      row++;
    }

    in.close();
  }
  return res;
}

const int LOOP_NUM = 30000;
const double T_CTRL = 0.001;
const double T_S = 0.02; //sampling parameter
const int T_HOR = 60; //time horizon

int main(){

    FILE *gid = popen("gnuplot -persist", "w");
    FILE *fp = fopen("state.dat", "w");

    MatrixXd Qmat = readMat("Qmat.dat", 4, 4);
    MatrixXd Pmat = readMat("Pmat.dat", 4, 4);
    MatrixXd Rmat = readMat("Rmat.dat", 2, 2);

    VectorXd x, u;
    x = VectorXd::Zero(4);
    //x(0) = 1;
    x(3) = M_PI;
    x(1) = 1;
    u = VectorXd::Zero(2);
    
    SAC sac(Qmat, Pmat, Rmat);

    double t = 0;
    int sim_loop, control_time;
    sim_loop = control_time = 0;
    
    fprintf(gid, "plot '-' using 1:2 with lines\n");
    while(sim_loop < LOOP_NUM){
        if(control_time >= T_S/T_CTRL){
            sac.Optimize(sim_loop*T_S, x);
            control_time = 0;
            // cout << "u_A = " << sac.u_A(0) << " " << sac.u_A(1) << endl;
            // cout << "tau_A:" << sac.tau_A << " duration:" << sac.duration << endl; 
        }
        u = sac.Control(control_time*T_CTRL);
        x += sac.state_eq(sim_loop*T_CTRL, x, u)*T_CTRL;
        //cout << "u1 = " << u(0) << " u2 = " << u(1) << endl;
        fprintf(fp, "%lf, %lf, %lf, %lf, %lf\n",sim_loop*T_CTRL, x(0), x(1), x(2), x(3));
        fprintf(gid, "%lf, %lf\n",x(0), x(1));
        control_time++;
        sim_loop++;
    }
    fprintf(gid, "e\n"); 
    fflush(gid);
    fclose(fp);
    pclose(gid);
    return 0;
}