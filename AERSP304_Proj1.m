%% AERSP 304 Project 1
% Payton Glynn, Craig Stenstrom, Nicholas Giampetro
%% Lagrange Point 2
clc;
clear;
close all;

load('EM_L2-304P1');   

pnts = 1000;    % number of descretizations of time

t = linspace(0,T,pnts);         % Set up our time step

options = odeset('reltol',1e-12,'abstol',1e-12);
[t,x] = ode45(@(t,x) odefun(t,x,MU1), t, x0, options);

% Convert to Inertial Frame
for i = 1:pnts
XI(i) = cos(t(i))*x(i,1)-sin(t(i))*x(i,2);
YI(i) = sin(t(i))*x(i,1)+cos(t(i))*x(i,2);
end

% Perturbed Solution
options = odeset('reltol',1e-12,'abstol',1e-12);
[t,xP] = ode45(@(t,x) odefun(t,x,MU1), t, x0+perturbation' , options);
xP = x - xP ;
for i = 1:pnts
    departX(i,1) = sqrt(xP(i,1)^2+xP(i,2)^2);
    departX(i,2) = sqrt(xP(i,3)^2+xP(i,4)^2);
end

% linearized
 xlin = 1.15568;
 ylin = 0;
    
 p1 = ((MU1+xlin)^2+ylin^2)^(1/2);
 p2 = ((1-MU1-xlin)^2+ylin^2)^(1/2);
    
Uxx = 1 - ((1-MU1)/(p1^3)) + (3*(1-MU1)*(xlin+MU1)^2)/(p1^5) + (3*MU1*(-xlin-MU1+1)^2)/(p2^5) - (MU1)/(p2^3);
Uyy = 1 - (1-MU1)/(p1^3) + (3*(1-MU1)*ylin^2)/(p1^5) - (MU1)/(p2^3) + (3*MU1*ylin^2)/(p2^5);
Uxy = (3*(1-MU1)*(xlin+MU1)*ylin)/(p1^5) - (3*MU1*ylin*(-xlin-MU1+1))/(p2^5);
 
A = [0 0 1 0; 0 0 0 1; Uxx Uxy 0 2; Uxy Uyy -2 0];
[eigvec eigval] = eig(A);
const = eigvec\perturbation;

tlin=zeros(pnts);
Xlin=zeros(pnts,4);
for i = 1:pnts
    tlin(i) = t(i);
    Xlin(i,:) = (const(1)*eigvec(:,1)*exp(eigval(1,1)*t(i))) + (const(2)*eigvec(:,2)*exp(eigval(2,2)*t(i))) + (const(3)*eigvec(:,3)*exp(eigval(3,3)*t(i))) + (const(4)*eigvec(:,4)*exp(eigval(4,4)*t(i)));
end

% linear pos calc
linPos=sqrt(Xlin(:,1).^2+Xlin(:,2).^2);

% plots nominal for L2 in B frame
figure
plot(x(:,1),x(:,2));
title('yN(t) vs. xN(t), Nominal, Lagrange Point 2');
xlabel('xN(t)');
ylabel('yN(t)');
ax = gca ;
exportgraphics(ax,'L2_Bframe.jpg')

% plots L2 in inertial frame
figure
plot(XI,YI)
title('Lagrange Point 2 in Inertial Frame')
xlabel('XN(t)');
ylabel('YN(t)');
ax = gca ;
exportgraphics(ax,'L2_Nframe.jpg')

% plots the perturbed case delta x for both position and velocity for L4
figure
subplot(2,1,1)
plot (t, departX(:,1))
title('Lagrange Point 2 perturbation position')
xlabel('t')
ylabel('delta x')
subplot(2,1,2)
plot (t, departX(:,2))
title('Lagrange Point 2 perturbation velocity')
xlabel('t')
ylabel('delta x')
ax = gca ;
exportgraphics(ax,'L2_Perturbed.jpg')

% linear comparison plot
figure
plot(t,linPos);
hold on 
plot(t,departX(:,1));
hold off
title('Lagrange Point 2 perturbation vs linearized motion')
xlabel('t')
ylabel('Pos')
legend('linear','perturbed')
ax = gca ;
exportgraphics(ax,'L2_Linear.jpg')


%% Lagrange Point 4

load('EM_L4-304P1')

t = linspace(0,T,pnts);         % Set up our time step

options = odeset('reltol',1e-12,'abstol',1e-12);
[t,x] = ode45(@(t,x) odefun(t,x,MU1), t, x0, options);

% Convert to Inertial Frame
for i = 1:pnts
XI(i) = cos(t(i))*x(i,1)-sin(t(i))*x(i,2);
YI(i) = sin(t(i))*x(i,1)+cos(t(i))*x(i,2);
end

% Perturbed Solution
options = odeset('reltol',1e-12,'abstol',1e-12);
[t,xP] = ode45(@(t,x) odefun(t,x,MU1), t, x0'+perturbation' , options);
xP = x - xP ;
for i = 1:pnts
    departX(i,1) = sqrt(xP(i,1)^2+xP(i,2)^2);
    departX(i,2) = sqrt(xP(i,3)^2+xP(i,4)^2);
end


% Linearize


% plots nominal for L4 in B frame
figure
plot(x(:,1),x(:,2))
title('yN(t) vs. xN(t), Nominal, Lagrange Point 4');
xlabel('xN(t)');
ylabel('yN(t)');
ax = gca ;
exportgraphics(ax,'L4_Bframe.jpg')

% plots L4 in inertial frame
figure
plot(XI,YI)
title('Lagrange Point 4 in Inertial Frame')
xlabel('XN(t)');
ylabel('YN(t)');
ax = gca ;
exportgraphics(ax,'L4_Nframe.jpg')

% plots the perturbed case delta x for both position and velocity for L4
figure
subplot(2,1,1)
plot (t, departX(:,1))
title('Lagrange Point 4 pertubation position')
xlabel('t')
ylabel('delta x')
subplot(2,1,2)
plot (t, departX(:,2))
title('Lagrange Point 4 pertubation velocity')
xlabel('t')
ylabel('delta x')
ax = gca ;
exportgraphics(ax,'L4_Perturbed.jpg')



%% Functions

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
