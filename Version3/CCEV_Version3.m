%% Electric Dc Engine control for a cruise control implementation
clc
clear all
close all

%% PARAMETERS
global rk lambda0 %Road Parameter selector
global K CL %CL closed loop simulation
global J m g rho S CD CR Km Bm r Rm Nf Nr Lm N
J = 1;                          % [kg/m^3] inertial (total, wheeel+Axes)
m = 1200;                       % [kg] Tesla model 3 mass
g = 9.81;                       % [m/s^2] gravity acceleration
CD = 2.3;                       % [-] drag coeffcient
rho = 1.225;                    % [kg/m^3] air desnity
S = 2.4;                        % [m^2] vehicle cross section
CR = 0.015;                     % [-] rolling resistance coefficient
Km = 0.75;                      % [N*m/A] Motor Torque Constant
Bm = 0.1;                       % [N*m*s/rad] Internal friction coefficent of the Motor
r = 0.25;                       % [m] Radious of the wheel
Rm = 0.4;                       % [ohm] Internal motor electric resistance
Lm = 0.05;                      % [H] reference armature inductance
Nf = 5880;                      % [N] Front normal force
Nr = 5880;                      % [N] Rear Normal force
N = Nr + Nf;                    % [N] Total Normal force
lambda0 = (0.0409821+0.13)/2;   % [-] reference slip ratio
rk = 1;                         % kind of road ... select one of the following
% 1) Dry asphalt
% 2) Wet asphalt
% 3) Snow
% 4) Ice
% 5) Dry Cobblestone
% 6) Wet cobblestone
% c = {'--b','--m','--c','--k','--g','--y','--r'};
global Thetaroad
Thetaroad = RoadCoefficents(rk);

%% Time span
global t0
t0 = 0; % [s] initial time
tf = 60; % [s] final time
dt = 0.1; % [s] time at which we want a solution 
tspan = (t0:dt:tf); % [s] time span :dt:

%% Linearisation Point

global xstar dstar fastar Imstar vstar omegastar Vmstar

%driving force
fastar = DrivingForceCCEV(rk,lambda0);

%disturbances
Lmstar = Lm;
CRstar = CR;
Nrstar = Nr;
Nfstar = Nf;
windstar = 0;
alphastar = 0;

dstar = [Lmstar, CRstar, Nfstar, Nrstar, windstar,alphastar ,Thetaroad]; %disturbance 
% in equilibrium condition supposed with no wind and no street inclination

%states

vstar = 130/3.6; % [m/s] linearisation speed
omegastar = vstar/(r*(1-lambda0));
Imstar = (1/Km)*(Bm*omegastar-(dstar(3)*dstar(2)-fastar)*r);

xstar = [vstar; omegastar; Imstar];

%control
Vmstar = Rm*Imstar + Km*omegastar;


%% Modal Analysis
% we consider the matrix A, B1, B2 and C calculated in the "PlantCCEV.m" file
global A B1 B2 C 
[A, B1, B2, C]= PlantCCEV2();
[W,D] = eig(A);
[V,Jordan] = jordan(A);

% Reachability
R = ctrb(A,B1);  % returns the reachability matrix [B AB A^2B .... A^(n-1)*B]

% Observability
O = obsv(A,C);   % returns the observability matrix [C; CA; CA^2 .... C*A^(n-1)]

res = TestStability(A); % Stability Check

D = 0;
StabDetCheck(A,B1,C,D); 

%% Observation Problem Setup

% Process noise
barQd11 = diag([0.01, 0.01, 0.01]);
% Mutual covariance
barQd12 = zeros(3,3);
% sensor noise
global std_tacho std_wheel std_amperometer
std_tacho = 1/3.6/3; % [m/s] standard deviation of the tachometer
std_wheel = 0.1/3; % [rad/s] standard deviation of the tone-wheel
std_amperometer = 0.02; % [A] standard deviation of the amperometer
% Measurement (sensor) noise
Rd = diag([std_tacho^2, std_wheel^2, std_amperometer^2]); 

%% Control Prolem Setup
% state penalization
v_max = 1/3.6; % [m/s] maximum allowable linear speed error
omega_max = 1; % [m] maximum allowable angular speed error
Im_max = 1; % [A] maximum allowable current error
Qp = inv(3*diag([v_max^2, omega_max^2, Im_max^2]));  
% input penalization
%u_max = 0.1*220; % [V] maximum allowable control action
u_max = 0.1;
Rp = inv(u_max^2);

%% LQG Solution

SYS = ss(A,B1,C,D); % declare the plant
QXU = blkdiag(Qp, Rp); % weigths associated to the control problem
QWV = [barQd11 barQd12; barQd12.' Rd]; % weigths associated to the observation problem
reg = lqg(SYS,QXU,QWV); % solution to the optimal control and observation problems
% matrices of the dynamic output feedback controller 
global H M KR
H = reg.A; 
M = reg.B;
KR = reg.C;

%% Open Loop Simulation 
CL = 0;  
% x0 = [zeros(6,1)]; % initial state 
x0 = [xstar; zeros(3,1)];
disp('STARTING OPEN LOOP SIMULATION')
[t,xol] = ode23(@LongDynamicsCCEV,tspan,x0); 
disp('ENDING OPEN LOOP SIMULATION')
% t = vector of time 
% xol = solution of the ode (open loop) (ode = ordinary differential equation)

%% Closed Loop Simulation
CL = 1;
disp('STARTING CLOSED LOOP SIMULATION')
[tcl,xcl] = ode23(@LongDynamicsCCEV,tspan,x0);
disp('ENDING CLOSED LOOP SIMULATION')

%% Plot 
figure(1)
plot(t,xol(:,1),'LineWidth',1.5)
xlabel('time [s]')
ylabel('v [m/s]')
title('Linear Speed')
grid on
hold on
plot(tcl,xcl(:,1),'LineWidth',1.5)
legend('Open Loop','Closed Loop')

figure(2)
plot(t,xol(:,2)*9.5493,'LineWidth',1.5)
xlabel('time [s]')
ylabel('\omega [rpm]')
title('Angular Speed')
grid on
hold on
plot(tcl,xcl(:,2),'LineWidth',1.5)
legend('Open Loop','Closed Loop')

figure(3)
plot(t,xol(:,3),'LineWidth',1.5)
xlabel('time [s]')
ylabel('Im [A]')
title('DC Engine Current')
grid on
hold on
plot(tcl,xcl(:,3),'LineWidth',1.5)
legend('Open Loop','Closed Loop')

%%
%
LT = [];
LC = [];
tick = 1;
%
n = 2;
assi = {'$t$ [s]','$\tilde{x}_1$ [m/s]','Linear speed error'};
legenda = {'Open Loop','Closed Loop'}; 
YMatrix = [(xol(:,1)-vstar).'; (xcl(:,1)-vstar).'];
createfigure(n, t, YMatrix, assi, legenda, tick, LT, LC)
%
n = 2;
assi = {'$t$ [s]','$\tilde{x}_2$ [rad/s]','Angular speed error'};
legenda = {'Open Loop','Closed Loop'}; 
YMatrix = [(xol(:,2)-omegastar).';(xcl(:,2)-omegastar).'];
createfigure(n, t, YMatrix, assi, legenda, tick, LT, LC)
%
n = 2;
assi = {'$t$ [s]','$\tilde{x}_3$ [A]','DC Current error'};
legenda = {'Open Loop','Closed Loop'}; 
YMatrix = [(xol(:,3)-Imstar).';(xcl(:,3)-Imstar).'];
createfigure(n, t, YMatrix, assi, legenda, tick, LT, LC)

% calculate the control action
u = zeros(length(tcl),1);
u0 = Rm*Imstar + Km*omegastar;
for i = 1:length(tcl)
 xi = xcl(i,4:6).';
 u(i) = u0 + KR*xi;
end
n = 1;
LT = {'-','--','--'};
assi = {'$t$ [s]','$\tilde{u}$ [V]','Control Voltage'};
legenda = {'LQG','$u_{max}$','$-u_{max}$'}; 
YMatrix = [(u-u0).'; u_max*ones(1,length(t)); -u_max*ones(1,length(t))];
createfigure(n, t, YMatrix, assi, legenda, tick, LT, LC)


