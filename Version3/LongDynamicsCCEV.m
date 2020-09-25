function dxdt = LongDynamicsCCEV(t,x)

%% PARAMETERS
global m g CD J rho S Bm Rm r Km dstar xstar rk lambda0 
Lm0 = dstar(1);
CR0 = dstar(2);
Nf0 = dstar(3);
Nr0 = dstar(4);
wind0 = dstar(5);
alpha0 = dstar(6);
theta1_0 = dstar(7);
theta2_0 = dstar(8);
theta3_0 = dstar(9);

v0 = xstar(1);
omega0 = xstar(2);
Im0 = xstar(3);

global CL H M KR

%% state assignment 

% x = col[v omega Im]   states vector
v = x(1);
omega = x(2);
Im = x(3);
xi = x(4:6);

if v > 36
    lambda = (x(2)*r-x(1))/(x(2)*r);
elseif v < 36
    lambda = (x(2)*r-x(1))/x(1);
end 

%% Disturbances and noises
[Lmn, CRn, Nfn, Nrn, w, alpha] = DisturbanceCCEV(t);  % disturbances
y = SensorsCCEV(t,x);  % measurements affected by noise
y0 = xstar;

%% CONTROL LAW
% u = Vmstar; % [V] DC Engine Armature Voltage
% u = ControlCCEV(t,x);  % control input
u0 = Rm*Im0 + Km*omega0;
if CL == 0 % open loop
    u = u0;
    dxidt = zeros(3,1);
elseif CL == 1 % closed loop
    dxidt = H*xi + M*(y-y0);  % observer
    u = u0 + KR*(x(1:3)-xstar);
end


%% DYNAMICS
dxdt = [(1/m)*(DrivingForceCCEV(rk,lambda)-m*g*sin(alpha)-0.5*rho*S*((v-w)^2)*CD-(Nr0+Nf0)*CR0);
        (1/J)*(-Bm*omega + Km*Im +(Nf0*CR0 - DrivingForceCCEV(rk,lambda))*r);
        (1/Lm0)*(-Km*omega - Rm*Im + u);
        dxidt];