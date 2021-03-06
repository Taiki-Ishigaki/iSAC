#include <iostream>
#include <fstream>
#include <string>
#include "sac.hpp"
#include "isac.hpp"

using namespace std;

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

static VectorXd state(double t, VectorXd x, VectorXd u){
    VectorXd dx(4);
    dx(0) =  u(0) * cos(x(3));
    dx(1) =  x(0) * sin(x(3));    
    dx(2) =  u(0);  
    dx(3) =  u(1);  
    return dx;
}

static MatrixXd control(double t, VectorXd x){
    MatrixXd h = MatrixXd::Zero(4,2);
    h(0, 0) =  cos(x(3)); 
    h(1, 0) =  sin(x(3)); 
    h(2, 0) = 1;
    h(3, 1) = 1;
    return h;
}

static MatrixXd dstate(double t, VectorXd x, VectorXd u){
    MatrixXd dfdx = MatrixXd::Zero(4,4);
    dfdx(0,3) =  -sin(x(3));     
    dfdx(1,3) =   cos(x(3));  
    return dfdx;
}

const int LOOP_NUM = 10000;
const double T_CTRL = 0.001;
const double T_S = 0.02; //sampling parameter
const int T_HOR = 60; //time horizon

int main(){

    FILE *gid = popen("gnuplot -persist", "w");
    FILE *fp = fopen("state.dat", "w");

    MatrixXd Qmat = readMat("Qmat.dat", 4, 4);
    MatrixXd Pmat = readMat("Pmat.dat", 4, 4);
    MatrixXd Rmat = readMat("Rmat.dat", 2, 2);

    VectorXd x, u, x_ref;
    x = VectorXd::Zero(4);
    //x(0) = 1;
    x(1) = 0;
    x(3) = 0;
    u = VectorXd::Zero(2);
    x_ref = VectorXd::Zero(4);
    x_ref(1) = -1;
    x_ref(0) = 1;
    x_ref(3) = 0;

    iSAC sac(4, 2);
    sac.Initialize(Qmat, Pmat, Rmat);
    double u_max[2] = { 1,  0.5};
    double u_min[2] = {-1, -0.5};
    sac.set_u_limit(u_max, u_min);

    sac.state_eq = state;
    sac.control_func = control;
    sac.dstate_eq = dstate;

    double t = 0;
    int sim_loop, control_time;
    sim_loop = control_time = 0;
    
    fprintf(gid, "plot '-' using 1:2 with lines\n");
    while(sim_loop < LOOP_NUM){
        if(control_time >= T_S/T_CTRL){
            sac.Optimize(sim_loop*T_S, x, x_ref);
            // x_ref(1) = sin(sim_loop*T_CTRL*2)/10;
            // x_ref(0) = sim_loop*T_CTRL/2;
            control_time = 0;
            cout << "u_A = " << sac.get_u_A() << endl;
            cout << "tau_A:" << sac.get_tau_A() << " duration:" << sac.get_duration() << endl; 
        }
        u = sac.Control(control_time*T_CTRL);
        x += sac.state_eq(sim_loop*T_CTRL, x, u)*T_CTRL;
        cout << "u = " << u(0)  << u(1)<< endl;
        cout << "x = " << x(0)  << x(1)<< endl;
        cout << "ref = " << x_ref(0) << x_ref(1) << endl;
        fprintf(fp, "%lf, %lf, %lf, %lf, %lf\n",sim_loop*T_CTRL, x(0), x(1), x(2), x(3));
        //fprintf(gid, "%lf, %lf\n",sim_loop*T_CTRL, x(0));
        fprintf(gid, "%lf, %lf\n",x(0), x(1));
        // fprintf(gid, "%lf, %lf\n",x_ref(0), x_ref(1));
        control_time++;
        sim_loop++;
    }
    fprintf(gid, "e\n"); 
    fflush(gid);
    fclose(fp);
    pclose(gid);
    return 0;
}