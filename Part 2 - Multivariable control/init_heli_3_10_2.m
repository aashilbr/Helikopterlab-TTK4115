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
%lambda_1 = -3+1i;
%lambda_2 = -3-1i;
%K_pp = (lambda_1 * lambda_2)/K_1;
%K_pd = -(lambda_1+lambda_2)/K_1;

%  load expoles1.mat;
%  load expoles2.mat;
%  load expoles3.mat;
%  load expoles4.mat;
 A = [ 0 1 0 0 0 ;
       0 0 0 0 0 ;
       0 0 0 0 0 ;
       1 0 0 0 0 ;
       0 0 1 0 0 ];
   
 A_poles = [0 1 0;
            0 0 0;
            0 0 0];
        
 B_poles = [0 0;
            0 K_1;
            K_2 0];
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
%test  
 Q = [125 0 0 0 0;
      0 80 0 0 0;
      0 0 150 0 0;
      0 0 0 1 0;
      0 0 0 0 1 ];
   
 R = [ 1.7 0 ;
      0 0.5 ];
  
 % p = [-1 -1+1i -1-1i];
% p = [-1 -2 -2];
% p = [-1 -2 -3];
  K = lqr(A,B,Q,R); 
  %K = place(A_poles,B_poles,p);
   disp(K);
  %F = inv(C*(inv(B*K-A))*B);
   F = [ K(1,1) K(1,3);
         K(2,1) K(2,3) ];
     
   %disp(F);
  %test
 
  load lqr1.mat;
  load lqr2.mat;
  load lqr6.mat;
  load lqr7.mat;
  load elv.mat;
  load pit.mat;
  load pit1.mat;
  load pit2.mat;
  load pit3.mat;
  load pit6.mat;
  load pit7.mat;
  load pit8.mat;
  load poles1.mat;
  load poles2.mat;
  load poles3.mat;
  load p_elevation1.mat;
  load p_elevation2.mat;
  load p_elevation3.mat;
  load p_pitch1.mat;
  load p_pitch2.mat;
  load p_pitch3.mat;
  load p_pitch4.mat;
  load lqri1.mat;
  load lqri2.mat;
  load lqri3.mat;
  load lqri4.mat;
  load lqri5.mat;
  load lqri6.mat;
  load lqri7.mat;
  load lqri8.mat;
  load lqri9.mat;
  load lqri10.mat;
  load lqri11.mat;
  load lqri12.mat;
  load lqri13.mat;
  load F_test1.mat;
  load F_test2.mat;
  load F_test3.mat;
  load evaluation.mat;
  load elv.mat;
  load pitch1.mat;
  load pitch2.mat;
  load pitch3.mat;
  
  figure(1);
  plot(pitch1(1,:), pitch1(9,:), 'LineWidth', 1.5, 'Color', 'k');
  hold on
  plot(pitch1(1,:), pitch1(4,:), 'LineWidth', 1.5, 'Color', 'r');
  plot(pitch2(1,:), pitch2(4,:), 'LineWidth', 1.5, 'Color', 'g');
  plot(pitch3(1,:), pitch3(4,:), 'LineWidth', 1.5, 'Color', 'b');
  title('Pitch angle');
  xlabel('time [s]');
 ylabel('Pitch angle [rad]');
 legend ('reference' , 'p1', 'p2', 'p3');
  % figure(1);
   %plot(elv3(1,:), elv3(8,:),'LineWidth', 1.5, 'Color', 'k');
 % hold on
 % plot(elv3(1,:), elv3(7,:),'LineWidth', 1.5, 'Color', 'b');
%   hold off
%   
%   figure(2);
%   plot(elv3(1,:), elv3(8,:),'LineWidth', 1.5, 'Color', 'k');
%   hold on
 %  plot(elv2(1,:), elv2(7,:),'LineWidth', 1.5, 'Color', 'r');
%   hold off
%   
%   figure(3);
%   plot(elv3(1,:), elv3(8,:),'LineWidth', 1.5, 'Color', 'k');
%   hold on
  % plot(elv1(1,:), elv1(7,:),'LineWidth', 1.5, 'Color', 'g');
%   hold off 
%   
%   figure(4);
%   plot(elv3(1,:), elv3(8,:),'LineWidth', 1.5, 'Color', 'k');
%   hold on
%   plot(elv4(1,:), elv4(7,:),'LineWidth', 1.5, 'Color', 'g');
%   hold off 
%   title('Elevation rate');
%  xlabel('time [s]');
% ylabel('Elevation rate [rad/s]');
%  legend ('reference' , 'LQR1', 'LQR2', 'LQR3');

%POLE PLACEMENT
 figure(1);
 %plot(p1(1,:), p1(8,:),'LineWidth', 1.5, 'Color', 'k')
%%plot(elv3(1,:), elv3(8,:),'LineWidth', 1.5, 'Color', 'k');
 %hold on
 %plot(p1(1,:), p1(7,:),'LineWidth', 1.5, 'Color', 'r')
 %plot(p2(1,:), p2(7,:),'LineWidth', 1.5, 'Color', 'b')
 %plot(p3(1,:), p3(7,:),'LineWidth', 1.5, 'Color', 'g')
% %plot(elv3(1,:), elv3(7,:),'LineWidth', 1.5, 'Color', 'b');
% hold off
% title('Elevation rate');
% xlabel('time [s]');
% ylabel('Elevation rate [rad/s]');
% legend ('reference' , 'p1', 'p2', 'p3');
% 
% 
%DENNE SKAL TESTES PÅ LAB, FREM TIL LQR
%  figure(2);
%  plot(pi1(1,:), pi1(9,:),'LineWidth', 1.5, 'Color', 'k')
%  hold on
%  plot(pi1(1,:), pi1(4,:),'LineWidth', 1.5, 'Color', 'r')
%  plot(pi2(1,:), pi2(4,:),'LineWidth', 1.5, 'Color', 'b')
%  plot(pi3(1,:), pi3(4,:),'LineWidth', 1.5, 'Color', 'g')
%  plot(pi4(1,:), pi4(4,:),'LineWidth', 1.5, 'Color', 'm')
%  hold off
%  title('Pitch angle');
%  xlabel('time [s]');
%  ylabel('Pitch angle [rad]');
%  legend ('reference' , 'p1', 'p2', 'p3');
%LQR WITHOUT INTEGRAL EFFECT
  %figure(1);
   %plot(poles1(1,:), poles1(8,:),'LineWidth', 1.5, 'Color', 'k');
   %hold on
 % plot(poles1(1,:), poles1(7,:), 'LineWidth', 1.5, 'Color', 'r');
   %plot(poles2(1,:), poles2(7,:), 'LineWidth', 1.5, 'Color', 'b');
   %plot(poles3(1,:), poles3(7,:), 'LineWidth', 1.5, 'Color', 'm');
 % hold off
  % plot(pit6(1,:), pit6(9,:),'LineWidth', 1.5, 'Color', 'k');
%  hold on
%plot(pit1(1,:), pit1(4,:),'LineWidth', 1.5, 'Color', 'r');
 %plot(pit2(1,:), pit2(4,:),'LineWidth', 1.5, 'Color', 'b');
%plot(pit3(1,:), pit3(4,:),'LineWidth', 1.5, 'Color', 'g');
%  %plot(pit4(1,:), pit4(4,:),'LineWidth', 1.5, 'Color', 'm'); dårlig
 %plot(pit5(1,:), pit5(4,:),'LineWidth', 1.5, 'Color', 'm');
 % plot(pit6(1,:), pit6(4,:),'LineWidth', 1.5, 'Color', 'm');
   %plot (pit7(1,:), pit7(4,:),'LineWidth', 1.5, 'Color', 'r');
   %plot (pit8(1,:), pit8(4,:),'LineWidth', 1.5, 'Color', 'b');
   %hold off
%  
%  %figure(2);
%   %plot(elv3(1,:), elv3(8,:),'LineWidth', 1.5, 'Color', 'k');
%  %hold on
%   %plot(elv5(1,:), elv5(7,:),'LineWidth', 1.5, 'Color', 'b');
%   %hold off

  % title('Pitch angle');
 %xlabel('time [s]');
 %ylabel('Pitch angle [rad]');
 %legend ('reference' , 'LQR1', 'LQR2', 'LQR3', 'LQR4');
  
%  %Plot elevation rate
%  figure(1);
%  %plot(Pol1(1,:), Pol1(8,:),'LineWidth', 1.5, 'Color', 'k');
%    plot(Pol4(1,:), Pol4(8,:),'LineWidth', 1.5, 'Color', 'k');
%  hold on
%  %plot(Pol1(1,:), Pol1(7,:),'LineWidth', 1.5, 'Color', 'b');
%  plot(Pol4(1,:), Pol4(7,:),'LineWidth', 1.5, 'Color', 'b');
% hold off

% 
% 
% %Plot pitch
% figure(2);
% plot(Pol4(1,:), Pol4(9,:),'LineWidth', 1.5, 'Color', 'k');
% hold on
% plot(Pol4(1,:), Pol4(4,:),'LineWidth', 1.5, 'Color', 'r');
% hold off
% title('Pitch angle');
% xlabel('time [s]');
% ylabel('Pitch angle [rad]');
% legend ('reference' , 'measured');
% %hold on
% %plot(Pol1(1,:), 

% %LQR WITH INTEGRAL EFFECT
% figure(1);
% plot(lqri1(1,:),lqri1(9,:),'LineWidth', 1.5, 'Color', 'k')
% hold on 
% plot(lqri1(1,:),lqri1(4,:),'LineWidth', 1.5, 'Color', 'r')
% plot(lqri2(1,:),lqri2(4,:),'LineWidth', 1.5, 'Color', 'b')
% plot(lqri3(1,:),lqri3(4,:),'LineWidth', 1.5, 'Color', 'g')
% plot(lqri4(1,:),lqri4(4,:),'LineWidth', 1.5, 'Color', 'm')
% %plot(lqri5(1,:),lqri5(4,:),'LineWidth', 1.5, 'Color', 'c')
% hold off
% figure(1);
% plot(lqri6(1,:),lqri6(9,:),'LineWidth', 1.5, 'Color', 'k')
% hold on
% plot(lqri6(1,:),lqri6(4,:),'LineWidth', 1.5, 'Color', 'r')
% plot(lqri7(1,:),lqri7(4,:),'LineWidth', 1.5, 'Color', 'y')
% plot(lqri8(1,:),lqri8(4,:),'LineWidth', 1.5, 'Color', 'g')
% plot(lqri9(1,:),lqri9(4,:),'LineWidth', 1.5, 'Color', 'm')
% %plot(lqri10(1,:),lqri10(4,:),'LineWidth', 1.5, 'Color', 'c')
% %plot(lqri11(1,:),lqri11(4,:),'LineWidth', 1.5, 'Color', 'y')
% plot(lqri12(1,:),lqri12(4,:),'LineWidth', 1.5, 'Color', 'c')
% plot(lqri13(1,:),lqri13(4,:),'LineWidth', 1.5, 'Color', 'b')
% 
% title('Pitch angle');
% xlabel('time [s]');
% ylabel('Pitch angle [rad]');
% legend ('reference' , 'LQR1', 'LQR2', 'LQR3', 'LQR4', 'LQR5', 'LQR6');
% hold off

% figure(1);
% plot(f2(1,:),f2(9,:),'LineWidth', 1.5, 'Color', 'k')
% hold on 
% plot(f2(1,:),f2(4,:),'LineWidth', 1.5, 'Color', 'r')
% plot(f3(1,:),f3(4,:),'LineWidth', 1.5, 'Color', 'b')
% plot(lqri12(1,:), lqri12(4,:),'LineWidth', 1.5, 'Color', 'g')
% hold off
% title('Pitch angle');
% xlabel('time [s]');
% ylabel('Pitch angle [rad]');
% legend ('reference' , 'F1', 'F2', 'F3');
% 
% figure(2);
% plot(evaluation(1,:),evaluation(9,:),'LineWidth', 1.5, 'Color', 'k')
% hold on 
% plot(evaluation(1,:),evaluation(4,:),'LineWidth', 1.5, 'Color', 'r')
% hold off
% figure(3);
% plot(evaluation(1,:),evaluation(8,:),'LineWidth', 1.5, 'Color', 'k')
% hold on 
% plot(evaluation(1,:),evaluation(5,:),'LineWidth', 1.5, 'Color', 'r')
% holf off