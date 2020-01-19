#ifndef _SAC_
#define _SAC_

#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

class SAC{
    private:

    public:
    /*parameter*/
    const double T_S = 0.03; //sampling parameter
    const int T_HOR = 10; //time horizon
    const double U_MAX = 500;
    const double U_MIN = -500;
    const double INF = 100000000.0;

    SAC(void);

    VectorXd state_eq(double t, VectorXd x, VectorXd u); //f
    MatrixXd control_func(double t, VectorXd x); //h
    MatrixXd dstate_eq(double t, VectorXd x, VectorXd u); //df

    double inc_cost(double t, VectorXd x);//incremental cost
    VectorXd dinc_cost(double t, VectorXd x);//differential incremental cost
    double end_cost(double t, VectorXd x); //end cost
    VectorXd dend_cost(double t, VectorXd x);//differential end cost

    VectorXd rho_eq(double t, VectorXd x_nom, VectorXd u_nom, VectorXd rho);

    double calc_J(double t, VectorXd x, VectorXd u);
    double calc_J_controlled(double t, MatrixXd x, MatrixXd u);

    void Optimize(double t, VectorXd x);

    VectorXd Control(double t);
    
    /*SAC parameter*/
    MatrixXd u_nom = VectorXd::Zero(2);//nominal control (often u_nom = 0)
    MatrixXd Q = MatrixXd::Identity(4,4);
    MatrixXd P = MatrixXd::Identity(4,4);
    MatrixXd R = MatrixXd::Identity(2,2);
    double alpha_d = -100;
    double default_duration = T_S/100;

    /*variable*/
    VectorXd u_A = u_nom;
    double tau_A;
    double duration;
};

#endif //_SAC_