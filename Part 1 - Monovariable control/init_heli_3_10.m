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
J_p = 2*m_p*(l_h)^2;
K_1 = K_f/(2*m_p*l_p);
lambda_1 = -3;
lambda_2 = -3;
K_pp = (lambda_1 * lambda_2)/K_1;
K_pd = -(lambda_1+lambda_2)/K_1;

%load conjugate.mat;
%load reell.mat;
%load imaginary.mat;
%load sammenfallende.mat;
%load sammenfallende1.mat;
%load sammenfallende2.mat;
%load sammenfallende3.mat;
%load conjugated1.mat;
%load conjugated2.mat;
%load conjugated3.mat;
load Reelle1.mat;
load Reelle2.mat;
load Reelle3.mat;
%load task1.mat
%plot(Conj(1,:),Conj(4,:), 'LineWidth', 1.5, 'Color', 'b');
%hold on
%plot(Ima(1,:), Ima(4,:), 'LineWidth', 1.5,'Color', 'r');
%hold on
%plot(Ree(1,:),Ree(4,:), 'LineWidth', 1.5,'Color', 'g');
%hold on%
%plot(Sam(1,:), Sam(4,:), 'LineWidth', 1.5,'Color', 'm');

%plot(Sam(1,:), Sam(8,:), 'LineWidth', 1.5, 'Color', 'k');
%plot(Conj1(1,:), Conj1(8,:),'LineWidth', 1.5, 'Color', 'k');
%plot(Sam1(1,:), Sam1(8,:), 'LineWidth', 1.5, 'Color', 'k');
%hold on
%plot(Sam1(1,:), Sam1(4,:), 'LineWidth', 1.5, 'Color', 'b');
%plot(Sam2(1,:), Sam2(4,:), 'LineWidth', 1.5, 'Color', 'r');
%plot(Sam3(1,:), Sam3(4,:), 'LineWidth', 1.5, 'Color', 'g');
plot(Ree1(1,:), Ree1(8,:), 'LineWidth', 1.5, 'Color', 'k');
hold on
plot(Ree1(1,:), Ree1(4,:), 'LineWidth', 1.5, 'Color', 'b');

plot(Ree2(1,:), Ree2(4,:), 'LineWidth', 1.5, 'Color', 'r');
plot(Ree3(1,:), Ree3(4,:), 'LineWidth', 1.5, 'Color', 'g');
hold off
title('Pitch angle');
xlabel('time [s]');
ylabel('pitch angle [rad]');
legend('reference', 'p1 = -1, p2 = -3', 'p1 = -9, p2 = -1', 'p1 = -20, p2 = -1');



%plot(Conj1(1,:), Conj1(4,:),'LineWidth', 1.5, 'Color', 'b');

%plot(Conj1(1,:), Conj1(8,:),'LineWidth', 1.5, 'Color', 'k');
%plot(Conj2(1,:), Conj2(4,:),'LineWidth', 1.5, 'Color', 'r');
%plot(Conj3(1,:), Conj3(4,:),'LineWidth', 1.5, 'Color', 'y');

%hold off
%title('Pitch angle');
%xlabel('time [s]');
%ylabel('pitch angle [rad]');
%legend ('reference' , 'p1 = -1, p2 = -3', 'p1 = -9, p2 = -1', 'p1 = -20, p2 = -1');

load Day1-Tuning.mat;

%plot(Tuning(1,:), Tuning(8,:),'LineWidth', 1.5, 'Color', 'k');
%hold on
%plot(Tuning(1,:), Tuning(4,:),'LineWidth', 1.5, 'Color', 'r');