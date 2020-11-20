% FOR HELICOPTER NR 3-10
% This file contains the initialization for the helicopter assignment in
% the course TTK4115. Run this file before you execute QuaRC_ -> Build 
% to build the file heli_q8.mdl.

% Oppdatert høsten 2006 av Jostein Bakkeheim
% Oppdatert høsten 2008 av Arnfinn Aas Eielsen
% Oppdatert høsten 2009 av Jonathan Ronen
% Updated fall 2010, Dominik Breu
% Updated fall 2013, Mark Haring
% Updated spring 2015, Mark Haring


%%%%%%%%%%% Calibration of the encoder and the hardware for the specific
%%%%%%%%%%% helicopter
Joystick_gain_x = 1;
Joystick_gain_y = -1;

%%%%%%%%
%DAY TWO
%%%%%%%%%%

%%%%%%%%%%% Physical constants
g = 9.81; % gravitational constant [m/s^2]
l_c = 0.46; % distance elevation axis to counterweight [m]
l_h = 0.66; % distance elevation axis to helicopter head [m]
l_p = 0.175; % distance pitch axis to motor [m]
m_c = 1.92; % Counterweight mass [kg]
m_p = 0.72; % Motor mass [kg]
V_s_0 = 7.8;
K_f = -(g*m_c*l_c-2*g*m_p*l_h)/(l_h*V_s_0);
J_p = 2*m_p*(l_h)^2;
J_e = m_c*(l_c)^2+2*m_p*(l_h)^2;
K_1 = K_f/(2*m_p*l_p);
K_2 = (l_h*K_f)/J_e;
lambda_1 = -3+1i;
lambda_2 = -3-1i;
K_pp = (lambda_1 * lambda_2)/K_1;
K_pd = -(lambda_1+lambda_2)/K_1;

%System matrices
 A = [ 0 1 0 0 0 ;
       0 0 0 0 0 ;
       0 0 0 0 0 ;
       1 0 0 0 0 ;
       0 0 1 0 0 ];
B = [ 0 0 ;
     0 K_1 ;
     K_2 0 ;
     0 0 ;
     0 0 ];
 
 C = [1 0 0 0 0;
      0 0 1 0 0];
  
  G = [ 0 0 ;
        0 0 ;
        0 0 ;
        1 0;
        0 1 ];
 
 A_poles = [0 1 0;
            0 0 0;
            0 0 0];
        
 B_poles = [0 0;
            0 K_1;
            K_2 0];

% LQR method
 Q = [125 0 0 0 0;
      0 80 0 0 0;
      0 0 150 0 0;
      0 0 0 1 0;
      0 0 0 0 1 ];
   
 R = [ 1.7 0 ;
      0 0.5 ];
 K = lqr(A,B,Q,R); 
 disp(K);
 %F = inv(C*(inv(B*K-A))*B);
  F = [ K(1,1) K(1,3);
         K(2,1) K(2,3) ];
%Pole placement  
% p = [-1 -1+1i -1-1i];
% p = [-1 -2 -2];
% p = [-1 -2 -3];
%K = place(A_poles,B_poles,p);
