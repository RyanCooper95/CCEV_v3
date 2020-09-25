%plant file for CCEV
function [A, B1, B2, C ] = PlantCCEV()
%partial derivatives evaluated from the "partialderivatives.m" script
%then replaced the simbolic values with the global ones
global Thetaroad Imstar vstar omegastar dstar Vmstar J m g rho S CD  Km Bm r Rm lambda0

a11 = (-1/m)*((dstar(3)/(omegastar*r))*(Thetaroad(1)*Thetaroad(2)*exp(-(lambda0*Thetaroad(2)))-Thetaroad(3))+rho*S*vstar*CD);
a12 = ((dstar(3)*vstar)/(m*r*omegastar^2))*((Thetaroad(1)*Thetaroad(2)*exp(-(lambda0*Thetaroad(2))))-Thetaroad(3));
a13 = 0;
a21 = (dstar(3)/(omegastar*J))*((Thetaroad(1)*Thetaroad(2)*exp(-(lambda0*Thetaroad(2))))-Thetaroad(3));
a22 = (-1/J)*(Bm+(vstar/omegastar^2)*(dstar(3)*Thetaroad(1)*Thetaroad(2)*exp(-(lambda0*Thetaroad(2)))-dstar(3)*Thetaroad(3)));
a23 = Km/J;
a31 = 0;
a32 = -Km/dstar(1);
a33 = -Rm/dstar(1);

b1_11 = 0;
b1_21 = 0;
b1_31 = 1/dstar(1);


b2_11 = 0;
b2_12 = -(dstar(3)+dstar(4))/m;
b2_13 = (Thetaroad(1)*(1-exp(-lambda0*Thetaroad(2)))-lambda0*Thetaroad(3)-dstar(2))/m;
b2_14 = -dstar(2)/m;
b2_15 = -rho*CD*S*(vstar-dstar(5))/m;
b2_16 = -g*cos(dstar(6));
b2_17 = (dstar(3)-dstar(4)*exp(-lambda0*Thetaroad(2)))/m;
b2_18 = (dstar(3)*Thetaroad(1)*Thetaroad(2)*exp(-lambda0*Thetaroad(2)))/m;
b2_19 = (-dstar(3)*lambda0)/m;
b2_21 = 0;
b2_22 = (dstar(3)*r)/J;
b2_23 = (dstar(2)-(Thetaroad(1)*(1-exp(-lambda0*Thetaroad(2)))-lambda0*Thetaroad(3)))*r/J;
b2_24 = 0;
b2_25 = 0;
b2_26 = 0;
b2_27 = (dstar(3)*exp(-lambda0*Thetaroad(2))-dstar(3))*r/J;
b2_28 = (-dstar(3)*Thetaroad(1)*Thetaroad(2)*exp(-lambda0*Thetaroad(2)))*r/J;
b2_29 = (dstar(3)*lambda0)*r/J;
b2_31 = (Km*omegastar - Vmstar + Imstar*Rm)/(dstar(1)^2);
b2_32 = 0;
b2_33 = 0;
b2_34 = 0;
b2_35 = 0;
b2_36 = 0;
b2_37 = 0;
b2_38 = 0;
b2_39 = 0;

c11 = 1;
c12 = 0;
c13 = 0;
c21 = 0;
c22 = 1;
c23 = 0;
c31 = 0;
c32 = 0;
c33 = 1;

A = [ a11 a12 a13
      a21 a22 a23
      a31 a32 a33];
  
B1 = [b1_11; b1_21; b1_31];

B2 = [b2_11 b2_12 b2_13 b2_14 b2_15 b2_16 b2_17 b2_18 b2_19
      b2_21 b2_22 b2_23 b2_24 b2_25 b2_26 b2_27 b2_28 b2_29
      b2_31 b2_32 b2_33 b2_34 b2_35 b2_36 b2_37 b2_38 b2_39];
  
C = [ c11 c12 c13
      c21 c22 c23
      c31 c32 c33];
  
  
  
end
