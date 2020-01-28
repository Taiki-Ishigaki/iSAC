#include "isac.hpp"

iSAC::iSAC(int x_dim, int u_dim){
    u_nom = VectorXd::Zero(u_dim);//nominal control (often u_nom = 0)
    u_A = u_nom;
    Q = MatrixXd::Identity(x_dim,x_dim);
    P = MatrixXd::Identity(x_dim,x_dim);
    R = MatrixXd::Identity(u_dim,u_dim);
    
    std::vector<VectorXd> u_s_(T_HOR, u_nom);
    std::vector<double> tau_s_(T_HOR, 0);
    std::vector<double> duration_s_(T_HOR, 0);
    u_s = u_s_;
    tau_s = tau_s_;
    duration_s = duration_s_;
}

void iSAC::Initialize(MatrixXd w1, MatrixXd w2, MatrixXd w3){
    Q  = w1;
    P  = w2;
    R  = w3;
}

void iSAC::set_u_limit(double max[], double min[]){
    u_max = max;
    u_min = min;
}

double iSAC::inc_cost(double t, VectorXd x, VectorXd x_ref){//incremental cost
    double cost = (x-x_ref).transpose()*Q*(x-x_ref);
    cost /= 2;
    return cost;
}
VectorXd iSAC::dinc_cost(double t, VectorXd x, VectorXd x_ref){//differential incremental cost
    VectorXd dcost = Q*(x-x_ref);
    return dcost;
}

double iSAC::end_cost(double t, VectorXd x, VectorXd x_ref){ //end cost
    double cost = (x-x_ref).transpose()*P*(x-x_ref);
    cost /= 2;
    return cost;
}

VectorXd iSAC::dend_cost(double t, VectorXd x, VectorXd x_ref){//differential end cost
    MatrixXd dcost = P*(x-x_ref);
    return dcost;
}

double iSAC::calc_J(double t, VectorXd x, VectorXd u, VectorXd x_ref){
    MatrixXd x_i[T_HOR+1];
    x_i[0] = x;
    double cost, integral_term, end_term;
    integral_term = inc_cost(t, x_i[0], x_ref);
    for(int loop_i = 1; loop_i <= T_HOR; loop_i++){
        x_i[loop_i] = x_i[loop_i-1] + state_eq(t + T_S*(loop_i-1), x_i[loop_i-1], u_def(t))*T_S;
        integral_term += inc_cost(t + T_S*loop_i, x_i[loop_i], x_ref);
    }
    end_term = end_cost(t + T_S*T_HOR, x_i[T_HOR], x_ref);
    cost = integral_term + end_term;
    return cost;
}

double iSAC::calc_J_controlled(double t, MatrixXd x, MatrixXd u, VectorXd x_ref){
    MatrixXd x_i[T_HOR+1], x_tmp;
    x_i[0] = x;
    double cost, integral_term, end_term;
    integral_term = inc_cost(t, x_i[0], x_ref);
    for(int loop_i = 1; loop_i <= T_HOR; loop_i++){
        if(((loop_i-1)*T_S < tau_A && tau_A < loop_i*T_S) &&
            ((loop_i-1)*T_S < tau_A+duration && tau_A+duration < loop_i*T_S)){
            x_tmp = x_i[loop_i-1] + state_eq(t + T_S*(loop_i-1), x_i[loop_i-1], u_def(T_S*(loop_i-1)))*(tau_A - (loop_i-1)*T_S);
            x_tmp = x_tmp + state_eq(tau_A, x_tmp, u_A)*duration;
            x_i[loop_i] = x_tmp + state_eq(tau_A+duration, x_tmp, u_def(tau_A+duration))*(loop_i*T_S-(tau_A+duration));
        }else if((loop_i-1)*T_S < tau_A && tau_A < loop_i*T_S){
            x_tmp = x_i[loop_i-1] + state_eq(t + T_S*(loop_i-1), x_i[loop_i-1], u_def(T_S*(loop_i-1)))*(tau_A - (loop_i-1)*T_S);
            x_i[loop_i] = x_tmp + state_eq(tau_A, x_tmp, u_A)*(loop_i*T_S - tau_A);
        }else if((loop_i-1)*T_S > tau_A && tau_A < loop_i*T_S){
            x_i[loop_i] = x_i[loop_i-1] + state_eq(t + T_S*(loop_i-1), x_i[loop_i-1], u_A)*T_S;
        }else if((loop_i-1)*T_S > tau_A && (loop_i-1)*T_S < tau_A+duration && tau_A+duration < loop_i*T_S){
            x_tmp = x_i[loop_i-1] + state_eq(t + T_S*(loop_i-1), x_i[loop_i-1], u_A)*(tau_A+duration - (loop_i-1)*T_S);
            x_i[loop_i] = x_tmp + state_eq(tau_A+duration, x_tmp, u_def(tau_A+duration))*(loop_i*T_S - tau_A+duration);
        }else{
            x_i[loop_i] = x_i[loop_i-1] + state_eq(t + T_S*(loop_i-1), x_i[loop_i-1], u_def(t))*T_S;
        } 
        integral_term += inc_cost(t + T_S*loop_i, x_i[loop_i], x_ref);
    }
    end_term = end_cost(t + T_S*T_HOR, x_i[T_HOR], x_ref);
    cost = integral_term + end_term;
    return cost;
}

void iSAC::Optimize(double t, VectorXd x, VectorXd x_ref){
    /*predict*/
    VectorXd x_nom[T_HOR+1];
    MatrixXd rho[T_HOR+1], jacob;
    x_nom[0] = x;
    for(int loop_i = 1; loop_i <= T_HOR; loop_i++){
        x_nom[loop_i] = x_nom[loop_i-1] + state_eq(t + T_S*loop_i, x_nom[loop_i-1], u_def(T_S*loop_i))*T_S;
    }
    rho[T_HOR] = dend_cost(t+T_S*T_HOR,x_nom[T_HOR],x_ref);
    MatrixXd E = MatrixXd::Identity(x.size(), x.size());
    for(int loop_i = T_HOR; loop_i > 0; loop_i--){
        jacob = dstate_eq(t + T_S*loop_i, x_nom[loop_i], u_def(T_S*loop_i));
        rho[loop_i-1] = (E -T_S*jacob.transpose()).inverse()*(rho[loop_i] + T_S*dinc_cost(t + T_S * (loop_i - 1), x_nom[loop_i-1], x_ref));
    }
    /* compute optimal action schedule u_s* */
    VectorXd u_opt_s[T_HOR+1];
    MatrixXd tmp_h_rho;
    for(int loop_i = 0; loop_i <= T_HOR; loop_i++){
        tmp_h_rho = control_func(t + T_S*loop_i, x_nom[loop_i]).transpose() * rho[loop_i];
        u_opt_s[loop_i] = u_def(T_S*loop_i) + (tmp_h_rho * tmp_h_rho.transpose() + R.transpose()).inverse() * (tmp_h_rho * tmp_h_rho.transpose()*u_def(T_S*loop_i) + tmp_h_rho*alpha_d);
    }
    /* Determine application time tau_A, input u_A */
    MatrixXd J_tau[T_HOR+1], J_taU_MIN;
    int min_index = 0;
    J_taU_MIN = rho[0].transpose()*(state_eq(t, x_nom[0], u_opt_s[0]) - state_eq(t, x_nom[0], u_def(0)));
    for(int loop_i = 1; loop_i <= T_HOR; loop_i++){
        J_tau[loop_i] = rho[loop_i].transpose()*(state_eq(t + T_S*loop_i, x_nom[loop_i], u_opt_s[loop_i]) - state_eq(t + T_S*loop_i, x_nom[loop_i], u_def(T_S*loop_i)));  
        if(J_taU_MIN(0,0) > J_tau[loop_i](0,0)){
            J_taU_MIN = J_tau[loop_i];
            min_index = loop_i;
        }
    };
    tau_A = T_S * min_index;
    u_A = u_opt_s[min_index];
    for(int loop_i = 0; loop_i < u_A.size(); loop_i++){
        if(u_A(loop_i) > u_max[loop_i]){
            u_A(loop_i) = u_max[loop_i];
        }else if(u_A(loop_i) < u_min[loop_i]){
            u_A(loop_i) = u_min[loop_i];
        }
    }
    /* Determine control duration lambda_A*/
    double J_new = INF;
    double J_init = calc_J(t, x, u_def(0), x_ref);
    double delta_J_min = -J_init;
    double dJ_prev;
    double omega = 1.2; //optimaize parametr
    double duration_prev;
    duration = duration_prev = default_duration;
    int k = 0;
    while(J_new -  J_init > delta_J_min && duration < T_S*T_HOR){
        duration = pow(omega, k) * duration;
        dJ_prev = J_new -  J_init;
        J_new = calc_J_controlled(t, x, u_def(0), x_ref);
        if(J_new -  J_init > dJ_prev){
            duration = duration_prev;
            break;
        }
        if(J_new -  J_init > 0){
            duration = default_duration;
            break;
        }
        duration_prev = duration;
        k += 1;
    }

    /*update stack*/
    u_s.push_back(u_A);
    u_s.erase(u_s.begin());
    tau_s.push_back(tau_A);
    tau_s.erase(tau_s.begin());
    duration_s.push_back(duration);
    duration_s.erase(duration_s.begin());
}

VectorXd iSAC::u_def(double t){
    VectorXd u = u_nom;
    for(int i = 0; i < T_HOR; i++){
        double late_time = T_S*(T_HOR-1-i);
        if(tau_s[i]-late_time <= t && t <= tau_s[i]-late_time + duration_s[i]-late_time){
            u = u_s[i];
        }
    }
    return u;
}

VectorXd iSAC::Control(double t){
    VectorXd u_out = u_def(t);
    if(tau_A <= t && t <= tau_A + duration){
        u_out = u_A;
    }
    return u_out;
}

VectorXd iSAC::get_u_A(void){
    return u_A;
}
double iSAC::get_tau_A(void){
    return tau_A;
}
double iSAC::get_duration(void){
    return duration;
}