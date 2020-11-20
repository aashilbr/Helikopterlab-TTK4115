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

L_1 = K_f *l_p;
L_2 = g*m_c*l_c-2*g*m_p*l_h;
L_3 = K_f * l_h;
L_4 = K_f*l_h;
J_p = 2*m_p*(l_p)^2;
J_e = m_c*(l_c)^2+2*m_p*(l_h)^2;
J_lambda = m_c*(l_c)^2+2*m_p*((l_h)^2+(l_p)^2);
K_1 = L_1/J_p;
K_2 = L_3/J_e;
K_3 = - (L_4*L_2)/(L_3*J_lambda);

PORT = 7;

%IMU OFFSET
IMU_offset = [0.018;
              -0.005;
              -0.09;
               0.01;
               0.004];

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
 
 A_kalman = [ 0 1 0 0 0 0 ;
       0 0 0 0 0 0;
       0 0 0 1 0 0;
       0 0 0 0 0 0;
       0 0 0 0 0 1;
       K_3 0 0 0 0 0];
B_kalman = [ 0 0 ;
     0 K_1 ;
     0 0;
     K_2 0 ;
     0 0 ;
     0 0 ];
 
 C_kalman = [1 0 0 0 0 0;
            0 1 0 0 0 0;
            0 0 1 0 0 0;
            0 0 0 1 0 0;
            0 0 0 0 0 1];
  
  G = [ 0 0 ;
        0 0 ;
        0 0 ;
        1 0;
        0 1 ];
%test  
%  Q = [80 0 0 0 0;
%       0 10 0 0 0;
%       0 0 50 0 0;
%       0 0 0 0.5 0;
%       0 0 0 0 0.5 ];
%    
%  R = [ 5 0 ;
%       0 10 ];
%   

Q = [125 0 0 0 0;
      0 80 0 0 0;
      0 0 150 0 0;
      0 0 0 1 0;
      0 0 0 0 1 ];
   
 R = [ 1.7 0 ;
      0 0.5 ];
  
  %p = [-1+3i -1 -1-3i];
  K = lqr(A,B,Q,R);
 % K = place(A,B,p);
   %disp(K);
  %F = inv(C*(inv(B*K-A))*B);
   F = [ K(1,1) K(1,3);
         K(2,1) K(2,3) ];
     
   %disp(F);
  % Q_d = [ 0.0001 0 0 0 0 0;
         %  0 0.0000003 0 0 0 0;
         %  0 0 0.00003 0 0 0;
          % 0 0 0 0.0000004 0 0;
          % 0 0 0 0 0 0;
          % 0 0 0 0 0 0.00003];
          
  Q_d = [0.1 0 0 0 0 0;
         0 0.0003 0 0 0 0;
         0 0 0.03 0 0 0;
         0 0 0 0.0004 0 0;
         0 0 0 0 0 0;
         0 0 0 0 0 0.03];

   load ground.mat;
   load linearization.mat;
   load switch.mat;
   load switchP.mat;
   load Qd1.mat;
   load Qd2.mat;
   load Qd3.mat;
   load Qd4.mat;
   %Task 1
   %Plot of helicopter laying still on the ground
 %   figure(1);
%    plot(ground(1,:), ground(6,:),'LineWidth', 1.5, 'Color', 'r');
%    title('Pitch');
%   xlabel('time [s]');
%  ylabel('Pitch angle [rad]');
%    figure(2);
%    plot(ground(1,:), ground(2,:),'LineWidth', 1.5, 'Color', 'r');
%    title('Pitch rate');
%    xlabel('time [s]');
%    ylabel('Pitch rate [rad/s]');
%    figure(3);
%    plot(ground(1,:), ground(3,:),'LineWidth', 1.5, 'Color', 'r');
%    title('Elevation');
%    xlabel('time [s]');
%    ylabel('Elevation [rad]');
%    figure(4);
%    plot(ground(1,:), ground(4,:), 'LineWidth', 1.5, 'Color', 'r');
%    title('Elevation rate');
%    xlabel('time [s]');
%    ylabel('Elevation rate [rad/s]');
%    figure(5);
%    plot(ground(1,:), ground(4,:), 'LineWidth', 1.5, 'Color', 'r');
%    title('Travel rate');
%    xlabel('time [s]');
%    ylabel('Travel rate [rad/s]');
%    
%    %Plots of linearization 
%     figure(6);
%     plot(linearization(1,1:10001), linearization(2,1:10001),'LineWidth', 1.5, 'Color', 'b');
%     title('Pitch ');
%     xlabel('time [s]');
%     ylabel('Pitch [rad]');
%    figure(7);
%    plot(linearization(1,:), linearization(3,:),'LineWidth', 1.5, 'Color', 'b');
%    title('Pitch rate');
%    xlabel('time [s]');
%    ylabel('Pitch rate [rad/s]');
%     figure(8);
%     plot(linearization(1,:), linearization(4,:),'LineWidth', 1.5, 'Color', 'b');
%     title('Elevation');
%     xlabel('time [s]');
%     ylabel('Elevation [rad]');
%     figure(9);
%     plot(linearization(1,:), linearization(5,:),'LineWidth', 1.5, 'Color', 'b');
%     title('Elevation rate');
%     xlabel('time [s]');
%     ylabel('Elevation rate [rad/s]');
%    figure(10);
%    plot(linearization(1,:), linearization(6,:),'LineWidth', 1.5, 'Color', 'b');
%    title('Travel rate');
%    xlabel('time [s]');
%    ylabel('Travel rate [rad/s]');
   %7501:17501
   
   T = transpose(ground(2:6,:));
   M1 = mean(T);
   C = cov(T)
   L = transpose(linearization(2:6,:));
   M2 = mean(L);
   lin = transpose(linearization(4,4966:10001));
   M3=mean(lin);
   R_d = cov(L)
   %Gjøre det samme for linearization
   %C = cov(transpose(ground(2:6));
   %disp(C);
   
   %Task 2
   T = 0.002;
   sys = ss(A_kalman,B_kalman,C_kalman,0);
   sys_d = c2d(sys, T);
   [A_d,B_d,C_d,~] = ssdata(sys_d);
  disp(A_d);
  disp(B_d);
  disp(C_d);
  disp(C);
  %plots of switch and x-hat
%   figure(2);
%   %plot(swi(1,:), swi(2,:),'LineWidth', 1.5, 'Color', 'k');
%   %hold on
%     plot(swi(1,:), swi(3,:),'LineWidth', 1.5, 'Color', 'b');
%     %hold off
%     title('Estimated state for pitch angle');
%     xlabel('time [s]');
%     ylabel('Pitch [rad]');
%     figure(3);
%     %plot(swi(1,:), swi(2,:),'LineWidth', 1.5, 'Color', 'k');
%     %hold on
%     plot(swi(1,:), swi(4,:),'LineWidth', 1.5, 'Color', 'b');
%     %hold off
%     title('Pitch rate');
%     xlabel('time [s]');
%     ylabel('Pitch rate [rad/s]');
%     legend('switch', 'pitch rate')
%     figure(4);
%     %plot(swi(1,:), swi(2,:),'LineWidth', 1.5, 'Color', 'k');
%     %hold on
%     plot(swi(1,:), swi(4,:),'LineWidth', 1.5, 'Color', 'k');
%     %hold off
%     title('Elevation');
%    xlabel('time [s]');
%     ylabel('Elevation [rad]');
%     figure(5);
%     %plot(swi(1,:), swi(2,:),'LineWidth', 1.5, 'Color', 'k');
%     %hold on
%     plot(swi(1,:), swi(5,:),'LineWidth', 1.5, 'Color', 'b');
%     %hold off
%     title('Elevation rate');
%    xlabel('time [s]');
%     ylabel('Elevation rate[rad/s]');
%     figure(6);
%     %plot(swi(1,:), swi(2,:),'LineWidth', 1.5, 'Color', 'k');
%     %hold on
%     plot(swi(1,:), swi(6,:),'LineWidth', 1.5, 'Color', 'b');
%     %hold off
%     title('Travel');
%    xlabel('time [s]');
%     ylabel('Travel [rad]');
%     figure(7);
%     %plot(swi(1,:), swi(2,:),'LineWidth', 1.5, 'Color', 'k');
%     %hold on
%     plot(swi(1,:), swi(7,:),'LineWidth', 1.5, 'Color', 'b');
%     %hold off
%     title('Travel rate');
%    xlabel('time [s]');
%     ylabel('Travel rate [rad/s]');
%     
%     %plots of switch and P-hat
%     figure(8);
%   %plot(swip(1,:), swip(2,:),'LineWidth', 1.5, 'Color', 'k');
%  % hold on
%     plot(swip(1,:), swip(3,:),'LineWidth', 1.5, 'Color', 'b');
%     %hold off
%     title('Estimated error convariance for pitch angle');
%     xlabel('time [s]');
%     ylabel('Pitch [rad]');
%     figure(9);
%     %plot(swip(1,:), swip(2,:),'LineWidth', 1.5, 'Color', 'k');
%     %hold on
%     plot(swip(1,:), swip(4,:),'LineWidth', 1.5, 'Color', 'b');
%     %hold off
%     title('Pitch rate');
%     xlabel('time [s]');
%     ylabel('Pitch rate [rad/s]');
%     legend('switch', 'pitch rate')
%     figure(10);
%     %plot(swip(1,:), swip(2,:),'LineWidth', 1.5, 'Color', 'k');
%     %hold on
%     plot(swip(1,:), swip(4,:),'LineWidth', 1.5, 'Color', 'b');
%     %hold off
%     title('Elevation');
%    xlabel('time [s]');
%     ylabel('Elevation [rad]');
%     figure(11);
%     %plot(swip(1,:), swip(2,:),'LineWidth', 1.5, 'Color', 'k');
%     %hold on
%     plot(swip(1,:), swip(5,:),'LineWidth', 1.5, 'Color', 'b');
%     %hold off
%     title('Elevation rate');
%    xlabel('time [s]');
%     ylabel('Elevation rate[rad/s]');
%     figure(12);
%     %plot(swip(1,:), swip(2,:),'LineWidth', 1.5, 'Color', 'k');
%     %hold on
%     plot(swip(1,:), swip(6,:),'LineWidth', 1.5, 'Color', 'b');
%     %hold off
%     title('Travel');
%    xlabel('time [s]');
%     ylabel('Travel [rad]');
%     figure(13);
%     %plot(swip(1,:), swip(2,:),'LineWidth', 1.5, 'Color', 'k');
%     %hold on
%     plot(swip(1,:), swip(7,:),'LineWidth', 1.5, 'Color', 'b');
%     %hold off
%     title('Travel rate');
%    xlabel('time [s]');
%     ylabel('Travel rate [rad/s]');

figure(1);
plot(qd1(1,:), qd1(2,:),'LineWidth', 1.5, 'Color', 'b');
hold on
plot(qd1(1,:), qd1(15,:),'LineWidth', 1.5, 'Color', 'r');

plot(qd1(1,:), qd1(9,:),'LineWidth', 1.5, 'Color', 'g');
hold off
title('Pitch angle');
xlabel('time [s]');
ylabel('Pitch [rad]');
legend('IMU', 'Encoder', 'Estimate');

figure(2);
plot(qd2(1,:), qd2(2,:),'LineWidth', 1.5, 'Color', 'b');
hold on
plot(qd2(1,:), qd2(15,:),'LineWidth', 1.5, 'Color', 'r');

plot(qd2(1,:), qd2(9,:),'LineWidth', 1.5, 'Color', 'g');
hold off
title('Pitch angle');
xlabel('time [s]');
ylabel('Pitch [rad]');
legend('IMU', 'Encoder', 'Estimate');

figure(3);
plot(qd3(1,:), qd3(2,:),'LineWidth', 1.5, 'Color', 'b');
hold on
plot(qd3(1,:), qd3(15,:),'LineWidth', 1.5, 'Color', 'r');

plot(qd3(1,:), qd3(9,:),'LineWidth', 1.5, 'Color', 'g');
hold off
title('Pitch angle');
xlabel('time [s]');
ylabel('Pitch [rad]');
legend('IMU', 'Estimate', 'Encoder');

figure(4);
plot(qd4(1,:), qd4(2,:),'LineWidth', 1.5, 'Color', 'b');
hold on
plot(qd4(1,:), qd4(15,:),'LineWidth', 1.5, 'Color', 'r');

plot(qd4(1,:), qd4(9,:),'LineWidth', 1.5, 'Color', 'g');
hold off
title('Pitch angle');
xlabel('time [s]');
ylabel('Pitch [rad]');
legend('IMU', 'Estimate', 'Encoder');

   
   