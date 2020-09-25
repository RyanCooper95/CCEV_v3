function u = ControlCCEV(t,x)

global Rm Imstar Km omegastar

%% nominal control (constant)
u0 = Rm*Imstar+Km*omegastar;
% if CL == 0
%     u = 0;
% else
%     u = -K*(X-X0); % LQR
% end

u = u0;