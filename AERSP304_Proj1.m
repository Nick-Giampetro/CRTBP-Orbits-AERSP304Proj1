%% AERSP 304 Project 1
% Payton Glynn, Craig Stenstrom, Nicholas Giampetro
%% Lagrange Point 2
clc;
clear;
close all;

load('EM_L2-304P1');   

pnts = 1000;    % number of descretizations of time

L2_t = linspace(0,T,pnts);         % Set up our time step

[L2_t, L2_x, L2_XI, L2_YI, L2_deltaXPos, L2_deltaXVel, L2_linPos, L2_linVel] = dataGen(L2_t, pnts, MU1, x0, perturbation);

% plots nominal for L2 in B frame
figure
plot(L2_x(:,1),L2_x(:,2));
title('yN(t) vs. xN(t), Nominal, Lagrange Point 2');
xlabel('xN(t)');
ylabel('yN(t)');
ax = gca ;
exportgraphics(ax,'L2_Bframe.jpg')

% plots L2 in inertial frame
figure
plot(L2_XI,L2_YI)
title('Lagrange Point 2 in Inertial Frame')
xlabel('XN(t)');
ylabel('YN(t)');
ax = gca ;
exportgraphics(ax,'L2_Nframe.jpg')

% plots L2 departure position and velocity vs time
figure
subplot(2,1,1)
plot(L2_t,L2_deltaXPos);      % Plots departure postion vs time
title('Lagrange Point 2, departure position vs time');
xlabel('time')
ylabel('departure position')
ax = gca ;
exportgraphics(ax,'L2_PerturbedPos.jpg')
subplot(2,1,2)
plot(L2_t,L2_deltaXVel);      % Plots departure velocity vs time
title('Lagrange 2, departure velocity vs time');
xlabel('time');
ylabel('departure velocity');
ax = gca ;
exportgraphics(ax,'L2_PerturbedVel.jpg')

% plots linearized solution and departure solution vs time
figure
subplot(2,1,1)
plot(L2_t,L2_linPos);
hold on 
plot(L2_t,L2_deltaXPos);
hold off
title('Lagrange Point 2 departure vs linearized Position')
xlabel('t')
ylabel('Pos')
legend('linear','departure')
ax = gca ;
exportgraphics(ax,'L2_LinearPos.jpg')
subplot(2,1,2)
plot(L2_t,L2_linVel);
hold on 
plot(L2_t,L2_deltaXVel);
hold off
title('Lagrange Point 2 departure vs linearized Velocity')
xlabel('t')
ylabel('Vel')
legend('linear','departure')
ax = gca ;
exportgraphics(ax,'L2_LinearVel.jpg')



%% Lagrange Point 4

load('EM_L4-304P1')

L4_t = linspace(0,T,pnts);         % Set up our time step

[L4_t, L4_x, L4_XI, L4_YI, L4_deltaXPos, L4_deltaXVel, L4_linPos, L4_linVel] = dataGen(L4_t, pnts, MU1, x0', perturbation);

% plots nominal for L4 in B frame
figure
plot(L4_x(:,1),L4_x(:,2))
title('yN(t) vs. xN(t), Nominal, Lagrange Point 4');
xlabel('xN(t)');
ylabel('yN(t)');
ax = gca ;
exportgraphics(ax,'L4_Bframe.jpg')

% plots L4 in inertial frame
figure
plot(L4_XI,L4_YI)
title('Lagrange Point 4 in Inertial Frame')
xlabel('XN(t)');
ylabel('YN(t)');
ax = gca ;
exportgraphics(ax,'L4_Nframe.jpg')

% plots departure postion and velocity vs time
figure
subplot(2,1,1)
plot(L4_t,L4_deltaXPos);      % Plots departure postion vs time
xlabel('time');
ylabel('departure position');
title('Lagrange point 4, departure position vs time');
ax = gca ;
exportgraphics(ax,'L4_PerturbedPos.jpg')
subplot(2,1,2)
plot(L4_t,L4_deltaXVel);      % Plots departure velocity vs time
xlabel('time');
ylabel('departure velocity');
title('Lagrange point 4, departure velocity vs time');
ax = gca ;
exportgraphics(ax,'L4_PerturbedVel.jpg')

% linear comparison plot L4
figure
subplot(2,1,1)
plot(L4_t,L4_linPos);
hold on 
plot(L4_t,L4_deltaXPos);
hold off
title('Lagrange Point 4 departure vs linearized Position')
xlabel('t')
ylabel('Pos')
legend('linear','perturbed')
ax = gca ;
exportgraphics(ax,'L4_LinearPos.jpg')
subplot(2,1,2)
plot(L4_t,L4_linVel);
hold on 
plot(L4_t,L4_deltaXVel);
hold off
title('Lagrange Point 4 departure vs linearized Velocity')
xlabel('t')
ylabel('Vel')
legend('linear','perturbed')
ax = gca ;
exportgraphics(ax,'L4_LinearVel.jpg')

%% Functions
function [t, x, XI, YI, deltaXPos, deltaXVel, linPos, linVel] = dataGen(t, pnts, MU1, x0, perturbation)
    options = odeset('reltol',1e-12,'abstol',1e-12);
    [t,x] = ode45(@(t,x) odefun(t,x,MU1), t, x0, options);

    % Convert to Inertial Frame
    for i = 1:pnts
    XI(i) = cos(t(i))*x(i,1)-sin(t(i))*x(i,2);
    YI(i) = sin(t(i))*x(i,1)+cos(t(i))*x(i,2);
    end

    % Perturbed Solution
    initP = x0+perturbation';
    options = odeset('reltol',1e-12,'abstol',1e-12);
    [t,xP] = ode45(@(t,x) odefun(t,x,MU1), t, initP, options);
    deltaX = xP - x;
    deltaXPos = sqrt((deltaX(:,1)).^2+(deltaX(:,2)).^2);
    deltaXVel = sqrt((deltaX(:,3)).^2+(deltaX(:,4)).^2);

    % linearized
    [linPos,linVel] = linearizer(initP(1), initP(2), t, MU1, pnts, perturbation) ;
end 


function    xdot = odefun(t,x,MU1)
    x1 = x(1);              
    x2 = x(2);            
    x3 = x(3);            
    x4 = x(4);              

    p1 = ((MU1+x1)^2+x2^2)^(1/2);
    p2 = ((1-MU1-x1)^2+x2^2)^(1/2);
    Ux = x1-(((1-MU1)*(x1+MU1))/(p1^3)) - MU1*(x1-1+MU1)/(p2^3);
    Uy = x2 - (((1-MU1)*x2)/(p1^3)) - (MU1*x2)/(p2^3);

    x1dot = x3;
    x2dot = x4;
    x3dot = 2*x4 + Ux;
    x4dot = -2*x3 + Uy;

    xdot = [x1dot; x2dot; x3dot; x4dot];
end

function [pos,vel] = linearizer(initX, initY, t, MU1, pnts, perturbation)
    xlin = initX;
    ylin = initY;
    
    p1 = ((MU1+xlin)^2+ylin^2)^(1/2);
    p2 = ((1-MU1-xlin)^2+ylin^2)^(1/2);
    
    Uxx = 1 - ((1-MU1)/(p1^3)) + (3*(1-MU1)*(xlin+MU1)^2)/(p1^5) + (3*MU1*(-xlin-MU1+1)^2)/(p2^5) - (MU1)/(p2^3);
    Uyy = 1 - (1-MU1)/(p1^3) + (3*(1-MU1)*ylin^2)/(p1^5) - (MU1)/(p2^3) + (3*MU1*ylin^2)/(p2^5);
    Uxy = (3*(1-MU1)*(xlin+MU1)*ylin)/(p1^5) - (3*MU1*ylin*(-xlin-MU1+1))/(p2^5);
 
    A = [0 0 1 0; 0 0 0 1; Uxx Uxy 0 2; Uxy Uyy -2 0];
    [eigvec, eigval] = eig(A);
    const = eigvec\perturbation;

    tlin=zeros(pnts);
    Xlin=zeros(pnts,4);
    
    for i = 1:pnts
        tlin(i) = t(i);
        Xlin(i,:) = (const(1)*eigvec(:,1)*exp(eigval(1,1)*t(i))) + (const(2)*eigvec(:,2)*exp(eigval(2,2)*t(i))) + (const(3)*eigvec(:,3)*exp(eigval(3,3)*t(i))) + (const(4)*eigvec(:,4)*exp(eigval(4,4)*t(i)));
    end

    % linear pos calc
    pos = sqrt(Xlin(:,1).^2+Xlin(:,2).^2);
    vel = sqrt(Xlin(:,3).^2+Xlin(:,4).^2);
end
