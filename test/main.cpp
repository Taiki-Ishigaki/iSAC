#include <iostream>
#include "sac.hpp"

using namespace std;

void Mat_size(MatrixXd m){
    cout << "rows = " << m.rows() << " cols == "<< m.cols() << endl;
}

void Vec_size(VectorXd v){
    cout << "size = " << v.size() << endl;
}

const int LOOP_NUM = 500;
const double T_CTRL = 0.001;
const double T_S = 0.1; //sampling parameter
const int T_HOR = 100; //time horizon

int main(){

    FILE *gid = popen("gnuplot -persist", "w");

    VectorXd x, u;
    x = VectorXd::Ones(4);
    x(0) = 100;
    u = VectorXd::Zero(2);
    
    SAC sac;

    double t = 0;
    int sim_loop, control_time;
    sim_loop = control_time = 0;
    
    fprintf(gid, "plot '-' with lines\n");
    while(sim_loop < LOOP_NUM){
        if(control_time > T_S/T_CTRL){
            sac.Optimize(sim_loop*T_S, x);
            control_time = 0;
        }
        u = sac.Control(control_time);
        x += sac.state_eq(sim_loop*T_CTRL, x, u)*T_CTRL;
        //cout << "u1 = " << u(0) << " u2 = " << u(1) << endl;
        //cout << "u_A = " << sac.u_A(0) << " " << sac.u_A(1) << endl;
        // Vec_size(sac.u_A);
        // cout << "tau_A:" << sac.tau_A << " duration:" << sac.duration << endl; 
        fprintf(gid, "%d, %lf\n",sim_loop, x(0));
        control_time++;
        sim_loop++;
    }
    fprintf(gid, "e\n"); 
    fflush(gid);

    pclose(gid);
    return 0;
}