%plant file for CCEV
function [A, B1, B2, C ] = PlantCCEV()
    global vstar omegastar Imstar dstar Vmstar
    global Km Rm Bm
    global m g CD rho S CR r J
    global Nf Nr N Cr lambda0 rk Thetaroad
    
    %% MATRIX A
    a11 = (1/m)*((-rho*S*CD*vstar)+((1/(r*omegastar))*(dstar(3)*Thetaroad(1)))-((1/(r*omegastar))*(dstar(3)*Thetaroad(1)*Thetaroad(2))*(exp((-Thetaroad(2)*omegastar*r+Thetaroad(2)*vstar)/(omegastar*r)))));
    a12 = (1/m)*((-vstar*dstar(3)*Thetaroad(2)/(r*omegastar^2))+(vstar*dstar(3)*Thetaroad(1)*Thetaroad(2)/(r*omegastar^2)*(exp((-Thetaroad(2)*omegastar*r+Thetaroad(2)*vstar)/(omegastar*r)))));
    a13 = 0;
   
    a21 = (1/J)*((((1/omegastar)*(dstar(3)*Thetaroad(1)*Thetaroad(2)))*(exp((-Thetaroad(2)*omegastar*r+Thetaroad(2)*vstar)/(omegastar*r))))-((1/omegastar)*(dstar(3)*Thetaroad(3))));
    a22 = (1/J)*((-Bm)-(((1/(omegastar^2))*(vstar*dstar(3)*Thetaroad(1)*Thetaroad(2)))*(exp((-Thetaroad(2)*omegastar*r+Thetaroad(2)*vstar)/(omegastar*r))))-((1/(omegastar^2))*(dstar(3)*Thetaroad(3)*vstar)));
    a23= (Km/J);
    
    a31 = 0;
    a32= (-Km/dstar(1));
    a33= (-Rm/dstar(1));
    
    %% MATRIX B1  
     b_1 = 0;
     b_2 = 0;
     b_3 = (1/dstar(1));

    %% MATRIX B2
    b11 = 0;
    b12 = -(1/m)*(dstar(3)+dstar(4));
    b13 = -CR/m;
    b14 = -CR/m;
    b15 = (1/m)*(rho*S*CD*vstar);
    b16 = -g;
    b17 = (1/m)*(dstar(3)*(1-exp((-Thetaroad(2)*omegastar*r+Thetaroad(2)*vstar)/(omegastar*r))));
    b18 = (1/m)*(dstar(3)*Thetaroad(1)*((-Thetaroad(2)*omegastar*r+Thetaroad(2)*vstar)/(omegastar*r))*exp((-Thetaroad(2)*omegastar*r+Thetaroad(2)*vstar)/(omegastar*r)));
    b19 = (-1/m)*(dstar(3)*((-Thetaroad(2)*omegastar*r+Thetaroad(2)*vstar)/(omegastar*r)));

    b21 = 0;
    b22 = (1/J)*(r*dstar(3));
    b23 = 0;
    b24 = (r/J)*(CR-DrivingForceCCEV(rk, lambda0));
    b25 = 0;
    b26 = 0;
    b27 = (1/J)*((r*dstar(3))*(-1+exp((-Thetaroad(2)*omegastar*r+Thetaroad(2)*vstar)/(omegastar*r))));
    b28 = (1/J)*((-r*dstar(3))*Thetaroad(1)*exp((-Thetaroad(2)*omegastar*r+Thetaroad(2)*vstar)/(omegastar*r))*((-Thetaroad(2)*omegastar*r+Thetaroad(2)*vstar)/(omegastar*r)));
    b29 = (r*dstar(3)/J)*((-Thetaroad(2)*omegastar*r+Thetaroad(2)*vstar)/(omegastar*r));

    b31 =  (-1/((dstar(1))^2))*(Vmstar-Rm*Imstar-Km*omegastar);
    b32 = 0;
    b33 = 0;
    b34 = 0;
    b35 = 0;
    b36 = 0;
    b37 = 0;
    b38 = 0;
    b39 = 0;

    %% RETURN 
    A = [a11 a12 a13;
      a21 a22 a23;
      a31 a32 a33];
  
    B1 = [b_1; b_2; b_3];
    
   B2 = [b11 b12 b13 b14 b15 b16 b17 b18 b19;
      b21 b22 b23 b24 b25 b26 b27 b28 b29;
      b31 b32 b33 b34 b35 b36 b37 b38 b39];       
    
      n = length(A);
    C = eye(n);
end
