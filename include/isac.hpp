#ifndef _ISAC_
#define _ISAC_

#include <iostream>
#include <vector> 
#include <Eigen/Dense>

using namespace Eigen;

class iSAC{
    public:
    iSAC(int x_dim, int u_dim);

    void Initialize(MatrixXd w1,MatrixXd w2,MatrixXd w3);
    void set_u_limit(double max[], double min[]);
    
    void Optimize(double t, VectorXd x, VectorXd x_ref);
    VectorXd Control(double t);

    VectorXd (*state_eq)(double t, VectorXd x, VectorXd u); //f
    MatrixXd (*control_func)(double t, VectorXd x); //h
    MatrixXd (*dstate_eq)(double t, VectorXd x, VectorXd u); //df

    VectorXd get_u_A(void);
    double get_tau_A(void);
    double get_duration(void);

    private:
    int x_dimention;
    int u_dimention;

    /*parameter*/
    const double T_S = 0.02; //sampling parameter
    const int T_HOR = 100; //time horizon
    const double INF = 100000000.0;
    const double EPS = 1.0e-3;

    double inc_cost(double t, VectorXd x, VectorXd x_ref);//incremental cost
    VectorXd dinc_cost(double t, VectorXd x, VectorXd x_ref);//differential incremental cost
    double end_cost(double t, VectorXd x, VectorXd x_ref); //end cost
    VectorXd dend_cost(double t, VectorXd x, VectorXd x_ref);//differential end cost

    double calc_J(double t, VectorXd x, VectorXd u, VectorXd x_ref);
    double calc_J_controlled(double t, MatrixXd x, MatrixXd u, VectorXd x_ref);

    VectorXd u_def(double t);
    
    /*iSAC parameter*/
    MatrixXd u_nom = VectorXd::Zero(2);//nominal control (often u_nom = 0)
    MatrixXd Q = MatrixXd::Identity(4,4);
    MatrixXd P = MatrixXd::Identity(4,4);
    MatrixXd R = MatrixXd::Identity(2,2);
    double alpha_d = -5;
    double default_duration = T_S/100;

    /*limit*/
    double* u_max;
    double* u_min;

    /*variable*/
    VectorXd u_A;
    double tau_A;
    double duration;

    std::vector<VectorXd> u_s;
    std::vector<double> tau_s;
    std::vector<double> duration_s;
};

#endif //_iSAC_