% Half-Car Suspension System Problem %

% A half-car suspension system was modeled, and from its free-body
% diagrams, governing equations were determined. Rearranging these
% equations, the equations of motion of the system were developed. The set
% of equations of motion were expressed in matrix form, and now the matrix
% form of the set of equations of motion are to be solved using numerical
% methods.

% Written by Aki Fujinawa, McGill University
% Date: December 12th, 2020

clc
clear
close all

% Parameters %

mc = 1000; % Vehicle chassis mass, sprung, kg
m1 = 50; % Front wheel mass, unsprung, kg
m2 = 50; % Rear wheel mass, unsprung, kg
k1 = 2000; % Front wheel spring, N/m
k2 = 500; % Spring between front wheel and vehicle chassis, N/m
k3 = 2000; % Rear wheel spring, N/m
k4 = 500; % Spring between rear wheel and vehicle chassis, N/m
c1 = 100; % Damper between front wheel and vehicle chassis, Ns/m
c2 = 100; % Damper between rear wheel and vehicle chassic, Ns/m
G = 0; % Center of gravity of system, m
L1 = 1; % Distance from center of gravity of front wheel to G, m
L2 = 1.5; % Distance from center of gravity of rear wheel to G, m
Ic = 200; % Moment of inertia of vehicle chassis, m^4
F1 = 2; % External force
w = 0.3; % Angular frequency of external force

% Mass Matrix %
 
M = [mc 0  0  0;
     0  m1 0  0;
     0  0  m2 0;
     0  0  0  Ic];

% Damping Matrix %

C = [c1+c2         -c1     -c2    L1*c1-L2*c2;
     -c1           c1      0      -L1*c1;
     -c2           0       c2     L2*c2;
     L1*c1-L2*c2   -L1*c1  L2*c2  (L1.^2)*c1+(L2.^2)*c2];

% Stiffness Matrix %

K = [k2+k4        -k2     -k4    L1*k2-L2*k4;
     -k2          k1+k2   0      -L1*k2;
     -k4          0       k3+k4  L2*k4;
     L1*k2-L2*k4  -L1*k2  L2*k4  (L1.^2)*k2+(L2.^2)*k4];

% VARIABLE TRANSFORMATION METHOD %

% Finding the resonant frequency of the system %

[U,wres] = (eig(K,M)); % Extracting modal matrix and diagonal matrix of eigenvalues

Ut = transpose(U);

alpha = 2; % Constant multiplier for mass matrix
beta = 0.5; % Constant multiplier for stiffness matrix

Ctemp = alpha.*M + beta.*K; % Damping matrix as a linear combination of mass and stiffness matrices
Cp = Ut * Ctemp * U; % Damping matrix in modal space
DR = zeros(4,1);

for i = 1:4
    DR(i,1) = Cp(i,i)/(2*sqrt(wres(i,i))); % Damping ratios
end

% At this point, we have enough information to decouple the system of
% equations and solve the individual, single DOF ODEs. 

% 4th ORDER RUNGE KUTTA METHOD % 

dt = 0.01; % Time step
N = 4000; % N/100 gives the number of seconds
t = 0:dt:N; % t goes from 0 to N/100 seconds.
normalized = zeros(4,length(t)); % zeroed matrix of normalized coordinates
 
for i = 1:4 % we need 4 iterations of this method as we have 4 ODEs to solve
    
    y0 = 0; % initial condition, time 0, displacement 0
    v0 = 0; % initial condition, time 0, velocity 0
    
    f1 = @(t,y,v) v; % function of velocity
    f2 = @(t,y,v) -Cp(i,i)*v -(wres(i,i))*y + Ut(i,2)*k1*0.001*sin(0.75*t) + Ut(i,3)*k3*0.001*sin(0.75*(t-pi/2)); % function of acceleration
    h = dt; % setting time step to variable h 
    
    y = zeros(length(t),1); % initializing zero vector for x
    v = zeros(length(t),1); % initializing zero vector for v
    y(1) = y0; % using initial condition for x
    v(1) = v0; % using initial condition for v

    for j = 1:length(t)-1 % 4th order Runge Kutta iterations
        
        k1y = f1(t(j),     y(j),         v(j));
        k1v = f2(t(j),     y(j),         v(j));
        
        k2y = f1(t(j)+h/2, y(j)+h*k1y/2, v(j)+h*k1v/2);
        k2v = f2(t(j)+h/2, y(j)+h*k1y/2, v(j)+h*k1v/2);
        
        k3y = f1(t(j)+h/2, y(j)+h*k2y/2, v(j)+h*k2v/2);
        k3v = f2(t(j)+h/2, y(j)+h*k2y/2, v(j)+h*k2v/2);
        
        k4y = f1(t(j)+h,   y(j)+h*k3y,   v(j)+h*k3v);
        k4v = f2(t(j)+h,   y(j)+h*k3y,   v(j)+h*k3v);
        
        y(j+1) = y(j) + (h/6)*(k1y + 2*k2y + 2*k2y + k4y);
        v(j+1) = v(j) + (h/6)*(k1v + 2*k2v + 2*k3v + k4v);

    end
    
    for k = 1:N
        normalized(i,k) = y(k); % adding eta values to the corresponding row of the normalized coordinate matrix
    end
    
end

xgen = U*normalized; % reverting values back to general coordinates
xplot = zeros(1,length(t)); 

% PLOTTING VERTICAL DISPLACEMENT OVER TIME OF EACH NODE IN SYSTEM

figure

for i = 1:4

    for k = 1:(length(t))
        xplot(1,k) = xgen(i,k);
    end
    
    subplot(2,2,i)
    plot(t,xplot,'b') % plot solution y(t)
    grid on
    
    xlabel('Time (s)');
    ylabel('Displacement (m)');
            
    if i == 1
        title('Vehicle Chassis Displacement');
    end
    if i == 2
        title('Front Wheel Displacement');
    end
    if i == 3
        title('Rear Wheel Displacement');
    end
    if i == 4
        title('Vehicle Chassis Angular Displacement');
        ylabel('Angular Displacement (rad)');
    end

    xlim([0,40]);
end

% PLOTTING TIME RESPONSE OF TRANSFER FUNCTION

figure

for i = 1:4
    
    sys = tf([1],[1 2*DR(i)*sqrt(wres(i,i)) (wres(i,i))]);

    t = 0:0.01:20;
    u = Ut(i,2)*0.001*sin(0.75*t)+Ut(i,3)*0.001*sin(0.75*(t-pi/3));

    subplot(2,2,i)
    lsim(sys,u,t) 
    grid on
            
    if i == 1
        title('Vehicle Chassis Response');
    end
    if i == 2
        title('Front Wheel Response');
    end
    if i == 3
        title('Rear Wheel Response');
    end
    if i == 4
        title('Vehicle Chassis Angular Response');
    end
    
end 

% PLOTTING FREQUENCY RESPONSE OF TRANSFER FUNCTION %

sos1 = [1];
sos2 = [1, 2*DR(1)*sqrt(wres(1,1)), (wres(1,1))];
[h,w] = freqz(sos2,sos1,'whole',2001);
figure
plot(w/pi,20*log10(abs(h)))
ax = gca;
ax.YLim = [-20 20];
ax.XTick = 0:.5:2;
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
title('Frequency Response of Car Suspension System')


% PLOTTING TIME RESPONSE OF TRANSFER FUNCTION VS INPUT AT VARYING FREQUENCIES

freq = [pi/6;pi/3;pi/2;pi];

for i = 1:4
    
    sys = tf([1],[1 2*DR(1)*sqrt(wres(1,1)) (wres(1,1))]);

    t = 0:0.01:20;
    u = Ut(1,2)*0.001*sin(freq(i)*t)+Ut(1,3)*0.001*sin(freq(i)*(t-pi/3));
    
    subplot(2,2,i)
    lsim(sys,u,t,'b') 
    grid on
            
    if i == 1
        title('Periodic excitation w=pi/6');
    end
    if i == 2
        title('Periodic excitation w=pi/3');
    end
    if i == 3
        title('Periodic excitation w=pi/2');
    end
    if i == 4
        title('Periodic excitation w=pi');
    end
end

