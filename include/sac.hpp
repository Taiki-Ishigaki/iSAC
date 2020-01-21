#ifndef _SAC_
#define _SAC_

#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

class SAC{
    private:

    public:
    /*parameter*/
    const double T_S = 0.02; //sampling parameter
    const int T_HOR = 60; //time horizon
    const double U_MAX[2] = { 10,  4};
    const double U_MIN[2] = {-10, -4};
    const double INF = 100000000.0;
    const double EPS = 1.0e-3;

    SAC(MatrixXd w1,MatrixXd w2,MatrixXd w3);

    VectorXd state_eq(double t, VectorXd x, VectorXd u); //f
    MatrixXd control_func(double t, VectorXd x); //h
    MatrixXd dstate_eq(double t, VectorXd x, VectorXd u); //df

    double inc_cost(double t, VectorXd x, VectorXd x_ref);//incremental cost
    VectorXd dinc_cost(double t, VectorXd x, VectorXd x_ref);//differential incremental cost
    double end_cost(double t, VectorXd x, VectorXd x_ref); //end cost
    VectorXd dend_cost(double t, VectorXd x, VectorXd x_ref);//differential end cost

    double calc_J(double t, VectorXd x, VectorXd u, VectorXd x_ref);
    double calc_J_controlled(double t, MatrixXd x, MatrixXd u, VectorXd x_ref);

    void Optimize(double t, VectorXd x, VectorXd x_ref);

    VectorXd Control(double t);
    
    /*SAC parameter*/
    MatrixXd u_nom = VectorXd::Zero(2);//nominal control (often u_nom = 0)
    MatrixXd Q = MatrixXd::Identity(4,4);
    MatrixXd P = MatrixXd::Identity(4,4);
    MatrixXd R = MatrixXd::Identity(2,2);
    double alpha_d = -1000;
    double default_duration = T_S/100;

    /*variable*/
    VectorXd u_A = u_nom;
    double tau_A;
    double duration;
};

#endif //_SAC_