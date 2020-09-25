function [Lmn, CRn, Nfn, Nrn, w, alpha] = DisturbanceCCEV(t)
%giving the time vector it returns me the disturbance vector. Be aware, in
%this version from all the disturbances the only one supposed to be
%variable is the wind. 
%Set the variable type to 0 or 1 in order to obtain a pulse or a step
%function

global Lm CR Nf Nr

type = 2; % [-] 0 = pulse, 1 = step

switch type
    case 0
        %% PULSE DISTURBANCE
        td = 20; % [s] time at which a step wind appears
        dt = 1; % [s] pulse duration
        if abs(t - td) <= dt 
            w = -40/3.6; % [m/s] wind speed opposite to the vehicle direction (from nose to tail)
        else
            w = 0; % [m/s] wind speed
        end
    case 1
        %% STEP DISTURBANCE
        td = 20; % [s] time at which a step wind appears
        if t > td
            w = -40/3.6; % [m/s] wind speed opposite to the vehicle direction (from nose to tail)
        else
            w = 0; % [m/s] wind speed
        end
    otherwise
        w = 0; % [m/s] wind speed
end
alpha = 3; % street inclination [deg]
CRn = CR;
Lmn = Lm;
Nfn = Nf;
Nrn = Nr;
end