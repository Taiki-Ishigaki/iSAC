/**
   Tank Controller
   @author Shin'ichiro Nakaoka
*/

#include <cnoid/SimpleController>
#include <cnoid/SpotLight>
#include <cnoid/RangeCamera>
#include <cnoid/RangeSensor>
#include <cnoid/Joystick>

#include "sac.hpp"
#include "isac.hpp"

using namespace std;
using namespace cnoid;

double RtoRad(Matrix3 R){
    if(R(2,2) = 1){
        if(R(0,0) != 0 && R(0,1) != 0){
            return atan2(R(0,1), R(0,0));
        }else if(R(0,0) == 0){
            return M_PI/2 * ((R(1,0) > 0)?1:-1);    
        }else if(R(1,0) == 0){
            return M_PI * ((R(0,0) > 0)?1:-1);  
        }
    }
    return 0;
}

class tank_controller : public SimpleController
{
    SimpleControllerIO* io;
    bool usePseudoContinousTrackMode;
    Link::ActuationMode turretActuationMode;
    Link* trackL;
    Link* trackR;
    Link* turretJoint[2];
    Link* tank_chassis;
    double qref[2];
    double qprev[2];
    double dt;

    Position T_i;// = io->body()->rootLink()->position();
    Vector3 p_i, pre_p_i;// = T.translation();
    Matrix3 R_i, pre_R_i;// = T.rotation();
            
    iSAC *sac;
    
    double u_max[2] = { 1,  0.5};
    double u_min[2] = {-1, -0.5};

    VectorXd x, u, x_ref;
    double t = 0;
    int sim_loop, control_time;

    const int LOOP_NUM = 30000;
    const double T_CTRL = 0.001;
    const double T_S = 0.02; //sampling parameter
    const int T_HOR = 600; //time horizon

    struct DeviceInfo {
        DevicePtr device;
        int buttonId;
        bool prevButtonState;
        bool stateChanged;
        DeviceInfo(Device* device, int buttonId)
            : device(device),
              buttonId(buttonId),
              prevButtonState(false),
              stateChanged(false)
        { }
    };
    vector<DeviceInfo> devices;
    SpotLightPtr spotLight;
    
    Joystick joystick;

public:
    static VectorXd f_state(double t, VectorXd x, VectorXd u){
        VectorXd dx(4);
        dx(0) =  u(0) * cos(x(3));
        dx(1) =  x(0) * sin(x(3));    
        dx(2) =  u(0);  
        dx(3) =  u(1);  
        return dx;
    }

    static MatrixXd h_control(double t, VectorXd x){
        MatrixXd h = MatrixXd::Zero(4,2);
        h(0, 0) =  cos(x(3)); 
        h(1, 0) =  sin(x(3)); 
        h(2, 0) = 1;
        h(3, 1) = 1;
        return h;
    }

    static MatrixXd df_dstate(double t, VectorXd x, VectorXd u){
        MatrixXd dfdx = MatrixXd::Zero(4,4);
        dfdx(0,3) =  -sin(x(3));     
        dfdx(1,3) =   cos(x(3));  
        return dfdx;
    }
    
    virtual bool initialize(SimpleControllerIO* io) override
    {

        MatrixXd Q(4,4), P(4,4), R(2,2);
        Q  << 30, 0, 0, 0,
               0, 30, 0, 0,
               0, 0, 0, 0,
               0, 0, 0, 5;
        P  <<  20, 0, 0, 0,
                0, 20, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 10;
        R  << 1, 0,
              0, 0.5;

        x = VectorXd::Zero(4);
        x(0) = 0.0;
        x(1) = 0;
        x(3) = 0;
        u = VectorXd::Zero(2);
        x_ref = VectorXd::Zero(4);
        x_ref(1) = 2;
        x_ref(0) = 1;
        x_ref(3) = 0;

        T_i = io->body()->rootLink()->position();
        pre_p_i = T_i.translation();
        pre_R_i = T_i.rotation();
        
        sac = new iSAC(4, 2);
        sac->Initialize(Q, P, R);
        sac->set_u_limit(u_max, u_min);
        sac->state_eq = &f_state;
        sac->control_func = &h_control;
        sac->dstate_eq = &df_dstate;

        sim_loop = control_time = 0;

        this->io = io;
        ostream& os = io->os();
        Body* body = io->body();

        io->setLinkInput(io->body()->rootLink(), LINK_POSITION);

        usePseudoContinousTrackMode = true;
        turretActuationMode = Link::ActuationMode::JOINT_TORQUE;
        for(auto opt : io->options()){
            if(opt == "wheels"){
                usePseudoContinousTrackMode = false;
            }
            if(opt == "velocity"){
                turretActuationMode = Link::ActuationMode::JOINT_VELOCITY;
            }
        }

        if(usePseudoContinousTrackMode){
            trackL = body->link("TRACK_L");
            trackR = body->link("TRACK_R");

        } else {
            trackL = body->link("WHEEL_L0");
            trackR = body->link("WHEEL_R0");
        }

        tank_chassis = body->link("CHASSIS");

        if(!trackL || !trackR){
            os << "The tracks are not found." << endl;
            return false;
        }

        if(usePseudoContinousTrackMode){
            trackL->setActuationMode(Link::JOINT_SURFACE_VELOCITY);
            trackR->setActuationMode(Link::JOINT_SURFACE_VELOCITY);
        } else {
            trackL->setActuationMode(Link::JOINT_VELOCITY);
            trackR->setActuationMode(Link::JOINT_VELOCITY);
        }
        io->enableOutput(trackL);
        io->enableOutput(trackR);
        
        turretJoint[0] = body->link("TURRET_Y");
        turretJoint[1] = body->link("TURRET_P");
        for(int i=0; i < 2; ++i){
            Link* joint = turretJoint[i];
            if(!joint){
                os << "Turret joint " << i << " is not found." << endl;
                return false;
            }
            qref[i] = qprev[i] = joint->q();
            joint->setActuationMode(turretActuationMode);
            io->enableIO(joint);
        }

        dt = io->timeStep();

        devices = {
            { body->findDevice<SpotLight>("Light"),    Joystick::A_BUTTON },
            { body->findDevice<RangeCamera>("Kinect"), Joystick::B_BUTTON },
            { body->findDevice<Camera>("Theta"),       Joystick::X_BUTTON },
            { body->findDevice<RangeSensor>("VLP-16"), Joystick::Y_BUTTON }
        };
        spotLight = dynamic_pointer_cast<SpotLight>(devices[0].device);

        // Turn on all the devices
        for(auto& device : devices){
            device.device->on(true);
            device.device->notifyStateChange();
        }

        return true;
    }

    virtual bool control() override
    {
        joystick.readCurrentState();
        
        double pos[2];
        for(int i=0; i < 2; ++i){
            pos[i] = joystick.getPosition(
                i==0 ? Joystick::L_STICK_H_AXIS : Joystick::L_STICK_V_AXIS);
            if(fabs(pos[i]) < 0.2){
                pos[i] = 0.0;
            }
        }
        /*set the velocity of each tracks*/
        if(usePseudoContinousTrackMode){
            double k = 1.0;
            trackL->dq_target() = k * (-2.0 * pos[1] + pos[0]);
            trackR->dq_target() = k * (-2.0 * pos[1] - pos[0]);
            // trackL->q() += 0.01;//= k * (-2.0 * pos[1] + pos[0]);
            // trackR->q() += 0.01;//= k * (-2.0 * pos[1] - pos[0]);
        } else {
            double k = 4.0;
            trackL->dq_target() = k * (-pos[1] + pos[0]);
            trackR->dq_target() = k * (-pos[1] - pos[0]);
        }

        if(control_time >= T_S/dt){
            sac->Optimize(sim_loop*T_S, x, x_ref);
            cout << "u_A = " << sac->get_u_A();
            cout << "tau_A:" << sac->get_tau_A() << " duration:" << sac->get_duration() << endl; 
            // cout << "dq =  "  <<  u(0) - 1.25 * u(1) << " : " << u(0) + 1.25 * u(1)<< endl;
            // cout << "u =  "  << u(0) << " : " << u(1) << endl;
            // cout << "x =  "  << x(0) << " : " << x(1) << " : " << x(2) << " : " << x(3) << endl;
            control_time = 0;
        }
        u = sac->Control(control_time*dt);
        trackL->dq_target() = u(0) - 1.25 * u(1);
        trackR->dq_target() = u(0) + 1.25 * u(1);

        // cout << "dq =  "  << (x(2) + u(0)*dt) + 1.25 * (x(3) + u(1)*dt) << " : " << (x(2) + u(0)*dt) - 1.25 * (x(3) + u(1)*dt);
        cout << "u =  "  << u(0) << " : " << u(1);
        cout << "x =  "  << x(0) << " : " << x(1) << " : " << x(2) << " : " << x(3) << endl;
        control_time++;
        sim_loop++;

        T_i = io->body()->rootLink()->position();
        p_i = T_i.translation();
        R_i = T_i.rotation();

        double w;
        Vector3d v;
        v = R_i.transpose() * (p_i - pre_p_i) / dt;
        w = (RtoRad(R_i) - RtoRad(pre_R_i))/dt;

        // cout << "v" <<  v(0) << endl;
        // cout << "w" <<  w << " : " << RtoRad(R_i) <<  endl;

        
        // cout << "p =  "  << p_i(0) << " : " << p_i(1) << " : " << p_i(2) << endl;
        // cout << "vx:" << vx << "vy:" << vy << "w" << RtoRad(R_i) << "R" << R_i(0,0) << endl;
        // cout << "t =  " << sim_loop*dt<< endl;
        // cout << trackR->dq() << " " << trackL->dq() << endl;
        // cout << "dv =  " << (x(2) - v(0))/dt<< " : " << (x(3) - w)/dt << endl;
        x(0) = p_i(0);
        x(1) = p_i(1);
        x(2) += v(0)*dt;//( trackR->dq() + trackL->dq()) * 1.50 / 2; // (dq_R + dq_L) * 2piR /2
        x(3) = RtoRad(R_i);//( trackR->dq() - trackL->dq()) * 1.50 / (2 * 1.25); // (dq_R - dq_L) * 2piR / (2 + 2piL)

        pre_p_i = p_i;
        pre_R_i = R_i;

        static const double P = 200.0;
        static const double D = 50.0;

        for(int i=0; i < 2; ++i){
            Link* joint = turretJoint[i];
            double pos = joystick.getPosition(
                i==0 ? Joystick::R_STICK_H_AXIS : Joystick::R_STICK_V_AXIS);
            if(fabs(pos) < 0.15){
                pos = 0.0;
            }

            if(turretActuationMode == Link::JOINT_VELOCITY){
                joint->dq_target() = pos;

            } else if(turretActuationMode == Link::JOINT_TORQUE){
                double q = joint->q();
                double dq = (q - qprev[i]) / dt;
                double dqref = 0.0;
                double deltaq = 0.002 * pos;
                qref[i] += deltaq;
                dqref = deltaq / dt;
                joint->u() = P * (qref[i] - q) + D * (dqref - dq);
                qprev[i] = q;
            }
        }

        for(auto& info : devices){
            if(info.device){
                bool stateChanged = false;
                bool buttonState = joystick.getButtonState(info.buttonId);
                if(buttonState && !info.prevButtonState){
                    info.device->on(!info.device->on());
                    stateChanged = true;
                }
                auto spotLight = dynamic_pointer_cast<SpotLight>(info.device);
                // if(spotLight){
                //     if(joystick.getPosition(Joystick::R_TRIGGER_AXIS) > 0.1){
                //         spotLight->setBeamWidth(
                //             std::max(0.1f, spotLight->beamWidth() - 0.001f));
                //         stateChanged = true;
                //     } else if(joystick.getButtonState(Joystick::R_BUTTON)){
                //         spotLight->setBeamWidth(
                //             std::min(0.7854f, spotLight->beamWidth() + 0.001f));
                //         stateChanged = true;
                //     }
                // }
                info.prevButtonState = buttonState;
                if(stateChanged){
                    info.device->notifyStateChange();
                }
            }
        }

        return true;
    }
};

CNOID_IMPLEMENT_SIMPLE_CONTROLLER_FACTORY(tank_controller)