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


%%%%%%%%%%% Physical constants
g = 9.81; % gravitational constant [m/s^2]
l_c = 0.46; % distance elevation axis to counterweight [m]
l_h = 0.66; % distance elevation axis to helicopter head [m]
l_p = 0.175; % distance pitch axis to motor [m]
m_c = 1.92; % Counterweight mass [kg]
m_p = 0.72; % Motor mass [kg]
V_s_0 = 7.8;
K_f = -(g*m_c*l_c-2*g*m_p*l_h)/(l_h*V_s_0);

L_1 = K_f * l_p;
L_2 = g*m_c*l_c - 2*g*m_p*l_h;
L_3 = K_f * l_h;
L_4 = K_f *l_h;
J_p = 2*m_p*(l_p)^2;
J_e = m_c*(l_c)^2+2*m_p*(l_h)^2;
J_lambda = m_c*(l_c)^2+2*m_p*((l_h)^2+(l_p)^2);
K_1 = L_1/J_p;
K_2 = L_3/J_e;
K_3 = - (L_4*L_2)/(L_3*J_lambda);

J_p = 2*m_p*(l_p)^2

%Device manager
PORT = 7;

% Offset from the IMU
IMU_offset = [0.018;
              -0.005;
              -0.09;
              0.01;
              0.004];

% System matrices
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
A_e = [0 1 0 0 0;
       0 0 0 0 0;
       0 0 0 1 0;
       0 0 0 0 0;
       K_3 0 0 0 0];

B_e = [0 0;
       0 K_1;
       0 0;
       K_2 0;
       0 0];
       
C = [1 0 0 0 0;
    0 1 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1];

% LQR method
Q = [125 0 0 0 0;
     0 80 0 0 0;
     0 0 150 0 0;
     0 0 0 1 0;
     0 0 0 0 1];
 R = [1.7 0;
      0 2];
 K = lqr(A,B,Q,R);
 F = [ K(1,1) K(1,3);
         K(2,1) K(2,3) ];
 % Eigenvalues of the controller
 P = eig(A-B*K);

%Pole placement
  r = 10;
  theta = 10*pi/180;
  p=[-r -r*(cos(theta)+1i*sin(theta)) -r*(cos(-theta)+1i*(sin(-theta))) -r*(cos(theta/2)+1i*(sin(theta/2))) -r*(cos(-theta/2)+1i*(sin(-theta/2)))];
L = place(A_e',C',p)';

disp(L);
