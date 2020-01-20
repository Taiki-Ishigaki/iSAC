#include "sac.hpp"

SAC::SAC(MatrixXd w1,MatrixXd w2,MatrixXd w3){
    Q  = w1;
    P  = w2;
    R  = w3;
}

VectorXd SAC::state_eq(double t, VectorXd x, VectorXd u){
    VectorXd dx(4);
    dx(0) =  x(2) * cos(x(3));
    dx(1) =  x(2) * sin(x(3));    
    dx(2) =  u(0);  
    dx(3) =  u(1);  
    return dx;
}

MatrixXd SAC::control_func(double t, VectorXd x){
    MatrixXd h = MatrixXd::Zero(4,2);
    h(2, 0) =  1; 
    h(3, 1) =  1; 
    return h;
}

MatrixXd SAC::dstate_eq(double t, VectorXd x, VectorXd u){
    MatrixXd dfdx = MatrixXd::Zero(4,4);
    dfdx(0,2) =  cos(x(3));
    dfdx(0,3) =  x(2) * -sin(x(3));    
    dfdx(1,2) =  sin(x(3));  
    dfdx(1,3) =  x(2) * cos(x(3));  
    return dfdx;
}

double SAC::inc_cost(double t, VectorXd x){//incremental cost
    double cost = x.transpose()*Q*x;
    cost /= 2;
    return cost;
}
VectorXd SAC::dinc_cost(double t, VectorXd x){//differential incremental cost
    VectorXd dcost = Q*x;
    return dcost;
}

double SAC::end_cost(double t, VectorXd x){ //end cost
    double cost = x.transpose()*P*x;
    cost /= 2;
    return cost;
}

VectorXd SAC::dend_cost(double t, VectorXd x){//differential end cost
    MatrixXd dcost = P*x;
    return dcost;
}

VectorXd SAC::rho_eq(double t, VectorXd x_nom, VectorXd u_nom, VectorXd rho){
    MatrixXd drho(4,1);
    drho = -(dinc_cost(t, x_nom).transpose() - dstate_eq(t, x_nom, u_nom).transpose())*rho;
    return drho;
}

double SAC::calc_J(double t, VectorXd x, VectorXd u){
    MatrixXd x_i[T_HOR+1];
    x_i[0] = x;
    double cost, integral_term, end_term;
    integral_term = inc_cost(t, x_i[0]);
    for(int loop_i = 1; loop_i <= T_HOR; loop_i++){
        x_i[loop_i] = x_i[loop_i-1] + state_eq(t + T_S*(loop_i-1), x_i[loop_i-1], u)*T_S;
        integral_term += inc_cost(t + T_S*loop_i, x_i[loop_i]);
    }
    end_term = end_cost(t + T_S*T_HOR, x_i[T_HOR]);
    cost = integral_term + end_term;
    return cost;
}

double SAC::calc_J_controlled(double t, MatrixXd x, MatrixXd u){
    MatrixXd x_i[T_HOR+1], x_tmp;
    x_i[0] = x;
    double cost, integral_term, end_term;
    integral_term = inc_cost(t, x_i[0]);
    for(int loop_i = 1; loop_i <= T_HOR; loop_i++){
        if(((loop_i-1)*T_S < tau_A && tau_A < loop_i*T_S) &&
            ((loop_i-1)*T_S < tau_A+duration && tau_A+duration < loop_i*T_S)){
            x_tmp = x_i[loop_i-1] + state_eq(t + T_S*(loop_i-1), x_i[loop_i-1], u_nom)*(tau_A - (loop_i-1)*T_S);
            x_tmp = x_tmp + state_eq(tau_A, x_tmp, u_A)*duration;
            x_i[loop_i] = x_tmp + state_eq(tau_A+duration, x_tmp, u_nom)*(loop_i*T_S-(tau_A+duration));
        }else if((loop_i-1)*T_S < tau_A && tau_A < loop_i*T_S){
            x_tmp = x_i[loop_i-1] + state_eq(t + T_S*(loop_i-1), x_i[loop_i-1], u_nom)*(tau_A - (loop_i-1)*T_S);
            x_i[loop_i] = x_tmp + state_eq(tau_A, x_tmp, u_A)*(loop_i*T_S - tau_A);
        }else if(tau_A < (loop_i-1)*T_S && tau_A < loop_i*T_S){
            x_i[loop_i] = x_i[loop_i-1] + state_eq(t + T_S*(loop_i-1), x_i[loop_i-1], u_A)*T_S;
        }else if((loop_i-1)*T_S < tau_A+duration && tau_A+duration < loop_i*T_S){
             x_tmp = x_i[loop_i-1] + state_eq(t + T_S*(loop_i-1), x_i[loop_i-1], u_A)*(tau_A+duration - (loop_i-1)*T_S);
             x_i[loop_i] = x_tmp + state_eq(tau_A+duration, x_tmp, u_nom)*(loop_i*T_S - tau_A+duration);
        }else{
            x_i[loop_i] = x_i[loop_i-1] + state_eq(t + T_S*(loop_i-1), x_i[loop_i-1], u_nom)*T_S;
        } 
        integral_term += inc_cost(t + T_S*loop_i, x_i[loop_i]);
    }
    end_term = end_cost(t + T_S*T_HOR, x_i[T_HOR]);
    cost = integral_term + end_term;
    return cost;
}

void SAC::Optimize(double t, VectorXd x){
    /*predict*/
    VectorXd x_nom[T_HOR+1];
    MatrixXd rho[T_HOR+1], jacob;
    x_nom[0] = x;
    for(int loop_i = 1; loop_i <= T_HOR; loop_i++){
        x_nom[loop_i] = x_nom[loop_i-1] + state_eq(t + T_S*loop_i, x_nom[loop_i-1], u_nom)*T_S;
    }
    rho[T_HOR] = dend_cost(t+T_S*T_HOR,x_nom[T_HOR]);
    MatrixXd E = MatrixXd::Identity(x.size(), x.size());
    MatrixXd K;
    for(int loop_i = T_HOR; loop_i > 0; loop_i--){
        jacob = dstate_eq(t + T_S*loop_i, x_nom[loop_i], u_nom);
        rho[loop_i-1] = (E-T_S*jacob.transpose()).inverse()*(rho[loop_i] + T_S*dinc_cost(t + T_S * (loop_i - 1), x_nom[loop_i-1]));
        K = (E-T_S*jacob.transpose()).inverse();
    }
    /* compute optimal action schedule u_s* */
    VectorXd u_opt_s[T_HOR+1];
    MatrixXd tmp_h_rho;
    for(int loop_i = 0; loop_i <= T_HOR; loop_i++){
        tmp_h_rho = control_func(t + T_S*loop_i, x_nom[loop_i]).transpose() * rho[loop_i];
        u_opt_s[loop_i] = u_nom + ((tmp_h_rho * tmp_h_rho.transpose() + R.transpose()).inverse() * tmp_h_rho)*alpha_d;
    }
    /* Determine application time tau_A, input u_A */
    MatrixXd J_tau[T_HOR+1], J_taU_MIN;
    int min_index = 0;
    J_taU_MIN = rho[0].transpose()*(state_eq(t, x_nom[0], u_opt_s[0]) - state_eq(t, x_nom[0], u_nom));
    for(int loop_i = 1; loop_i <= T_HOR; loop_i++){
        J_tau[loop_i] = rho[loop_i].transpose()*(state_eq(t + T_S*loop_i, x_nom[loop_i], u_opt_s[loop_i]) - state_eq(t + T_S*loop_i, x_nom[loop_i], u_nom));  
        if(J_taU_MIN(0,0) > J_tau[loop_i](0,0)){
            J_taU_MIN = J_tau[loop_i];
            min_index = loop_i;
        }
    };
    //std::cout<< u_A << std::endl;
    tau_A = T_S * min_index;
    u_A = u_opt_s[min_index];
    for(int loop_i = 0; loop_i < u_A.size(); loop_i++){
        if(u_A(loop_i) > U_MAX[loop_i]){
            u_A(loop_i) = U_MAX[loop_i];
        }else if(u_A(loop_i) < U_MIN[loop_i]){
            u_A(loop_i) = U_MIN[loop_i];
        }
    }
    /* Determine control duration lambda_A*/
    double J_new = INF;
    double J_init = calc_J(t, x, u_nom);
    double delta_J_min = -J_init*0.1;
    double dJ_prev;
    double omega = 1.2; //optimaize parametr
    double duration_prev;
    duration = duration_prev = default_duration;
    int k = 0;
    while(J_new -  J_init > delta_J_min && duration < T_S*T_HOR){
        duration = pow(omega, k) * duration;
        dJ_prev = J_new -  J_init;
        J_new = calc_J_controlled(t, x, u_nom);
        if(J_new -  J_init > dJ_prev){
            duration = duration_prev;
            break;
        }
        if(J_new -  J_init > 0){
            duration = default_duration;
        }
        duration_prev = duration;
        k += 1;
    }
}

VectorXd SAC::Control(double t){
    if(t < tau_A){
        return u_nom;
    }else if(tau_A <= t && t <= tau_A + duration){
        return u_A;
    }else if(tau_A + duration < t){
        return u_nom;
    }else{
        return u_nom; 
    }
}