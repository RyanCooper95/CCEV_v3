function y = SensorsCCEV(t,x)

global std_tacho std_wheel std_amperometer
% nu = noises
nu_t = 0*std_tacho*randn(1);
nu_w = 0*std_wheel*randn(1);
nu_a = 0*std_amperometer*randn(1);
y = zeros(3,1);
y(1) = [1 0 0 0 0 0]*x+nu_t; 
y(2) = [0 1 0 0 0 0]*x+nu_w;
y(3) = [0 0 1 0 0 0]*x+nu_a;

end