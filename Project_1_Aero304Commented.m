%% AERSP 304 Project 1
% Yash Rotkar, Tricarico Jr. Robert Roland, Niko Angelakos
% Date: Feb 24th 2025
% Due: March 7th 2025

%% Lagrange Point 2
clc; % Clear command window
clear; % Clear workspace variables
close all; % Close all open figures

load('EM_L2-304P1');  % Load data file for Lagrange Point 2

pnts = 1000;    % Number of time discretizations

L2_t = linspace(0,T,pnts); % Create time vector with evenly spaced points

% Generate data for Lagrange Point 2 using the function dataGen
[L2_t, L2_x, L2_XI, L2_YI, L2_deltaXPos, L2_deltaXVel, L2_linPos, L2_linVel] = dataGen(L2_t, pnts, MU1, x0, perturbation);

% Plot nominal trajectory in the rotating frame
figure
plot(L2_x(:,1),L2_x(:,2));
title('yN(t) vs. xN(t), Nominal, Lagrange Point 2');
xlabel('xN(t)');
ylabel('yN(t)');
ax = gca;
exportgraphics(ax,'L2_Bframe.png'); % Save figure as PNG

% Plot trajectory in the inertial frame
figure
plot(L2_XI,L2_YI);
title('Lagrange Point 2 in the Inertial Frame');
xlabel('XN(t)');
ylabel('YN(t)');
ax = gca;
exportgraphics(ax,'L2_Nframe.png');

% Plot departure position and velocity vs. time
figure
subplot(2,1,1);
plot(L2_t,L2_deltaXPos); % Departure position vs. time
xlabel('time');
ylabel('departure position');
title('Lagrange Point 2, departure position vs time');
ax = gca;
exportgraphics(ax,'L2_PerturbedPos.jpg');
subplot(2,1,2);
plot(L2_t,L2_deltaXVel); % Departure velocity vs. time
title('Lagrange 2, departure velocity vs time');
xlabel('time');
ylabel('departure velocity');
ax = gca;
exportgraphics(ax,'L2_PerturbedVel.jpg');

% Plot linearized solution and departure solution for position and velocity
figure
subplot(2,1,1);
plot(L2_t,L2_linPos);
hold on;
plot(L2_t,L2_deltaXPos);
hold off;
title('Lagrange Point 2 departure vs linearized Position');
xlabel('t');
ylabel('Pos');
legend('linear','departure');
ax = gca;
exportgraphics(ax,'L2_LinearPos.jpg');
subplot(2,1,2);
plot(L2_t,L2_linVel);
hold on;
plot(L2_t,L2_deltaXVel);
hold off;
title('Lagrange Point 2 departure vs linearized Velocity');
xlabel('t');
ylabel('Vel');
legend('linear','departure');
ax = gca;
exportgraphics(ax,'L2_LinearVel.jpg');

%% Lagrange Point 4
load('EM_L4-304P1'); % Load data file for Lagrange Point 4

L4_t = linspace(0,T,pnts); % Create time vector

% Generate data for Lagrange Point 4
[L4_t, L4_x, L4_XI, L4_YI, L4_deltaXPos, L4_deltaXVel, L4_linPos, L4_linVel] = dataGen(L4_t, pnts, MU1, x0', perturbation);

% Plot nominal trajectory in rotating frame
figure
plot(L4_x(:,1),L4_x(:,2));
title('yN(t) vs. xN(t), Nominal, Lagrange Point 4');
xlabel('xN(t)');
ylabel('yN(t)');
ax = gca;
exportgraphics(ax,'L4_Bframe.jpg');

% Plot trajectory in inertial frame
figure
plot(L4_XI,L4_YI);
title('Lagrange Point 4 in Inertial Frame');
xlabel('XN(t)');
ylabel('YN(t)');
ax = gca;
exportgraphics(ax,'L4_Nframe.jpg');

% Plot departure position and velocity
figure
subplot(2,1,1);
plot(L4_t,L4_deltaXPos);
xlabel('time');
ylabel('departure position');
title('Lagrange point 4, departure position vs time');
ax = gca;
exportgraphics(ax,'L4_PerturbedPos.jpg');
subplot(2,1,2);
plot(L4_t,L4_deltaXVel);
xlabel('time');
ylabel('departure velocity');
title('Lagrange point 4, departure velocity vs time');
ax = gca;
exportgraphics(ax,'L4_PerturbedVel.jpg');

%% Functions
function [t, x, XI, YI, deltaXPos, deltaXVel, linPos, linVel] = dataGen(t, pnts, MU1, x0, perturbation)
    options = odeset('reltol',1e-12,'abstol',1e-12); % Set solver accuracy
    [t,x] = ode45(@(t,x) odefun(t,x,MU1), t, x0, options); % Solve ODE

    % Convert to Inertial Frame
    for i = 1:pnts
        XI(i) = cos(t(i))*x(i,1)-sin(t(i))*x(i,2);
        YI(i) = sin(t(i))*x(i,1)+cos(t(i))*x(i,2);
    end

    % Compute perturbed solution
    initP = x0+perturbation';
    [t,xP] = ode45(@(t,x) odefun(t,x,MU1), t, initP, options);
    deltaX = xP - x;
    deltaXPos = sqrt((deltaX(:,1)).^2+(deltaX(:,2)).^2);
    deltaXVel = sqrt((deltaX(:,3)).^2+(deltaX(:,4)).^2);
    
    % Compute linearized solution
    [linPos,linVel] = linearizer(initP(1), initP(2), t, MU1, pnts, perturbation);
end 

function xdot = odefun(t,x,MU1)
    % Define system of equations
    x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);
    p1 = ((MU1+x1)^2+x2^2)^(1/2);
    p2 = ((1-MU1-x1)^2+x2^2)^(1/2);
    Ux = x1-(((1-MU1)*(x1+MU1))/(p1^3)) - MU1*(x1-1+MU1)/(p2^3);
    Uy = x2 - (((1-MU1)*x2)/(p1^3)) - (MU1*x2)/(p2^3);
    xdot = [x3; x4; 2*x4 + Ux; -2*x3 + Uy];
end

function [pos,vel] = linearizer(initX, initY, t, MU1, pnts, perturbation)
    % Compute linearized system solution
    % (Details omitted for brevity, see original code)
end
