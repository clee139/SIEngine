function Team05_Project2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project #2
% Payton Downey, Charles Lee, Kade Pizzuto, Nicholas Tate, Jaylin Trice, Ethan Winchester
% ME 2543--Simulations Methods
% Spring 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; % Clear the variable list and the command window.
% Problem Number: 1
gamma = 1.4; % Specific Heat Ratio of the Gas
R = 287; % Gas Constant (J/kg/K)
r = 8.4; % Compression ratio of the engine
l = 0.12; % Connecting rod length (m)
S = 0.08; % Piston Stroke (m)
V0 = 5.0e-5; % Cylinder clearance volume (m^3)
b = 0.09; % Cylinder bore (diameter) (m)
T1 = 300; % Inital Temperature (K)
theta_s = (3*pi)/2; % Crank angle at the start of heat addition
theta_b = pi; % Angular interval over which heat is added
q_in = 2.8e6; % Total mass-specific heat added to the gas due to combustion (J/kg)
Tw = 300; % Cylinder wall temperature (K)
P1 = 1.013e5; % Pressure (Pa) 
W_in = 0; % Work in (J)
C = 0; % Empirical proportionality constant (s^-1)
h = 0; % Convective heat transfer parameter (W/m^2K)
epsilon = S/(2*l);
M0 = (P1*V0)/(T1*R);
w = 50; % Engine speed (rad/s); Irreelevant for Problem 1
y0 = [P1;W_in]; % Initial Guess 
theta_domain = [pi (3*pi)]; % Range for Theta Values
% Solve for Pressure and Work In using ode45
[theta,sol] = ode45(@myeqs,theta_domain,y0);
% Solve for Volume 
[Vol1,dVol1] = myv(theta);
% Solve for Temperature 
Temp1 = myt(Vol1,mym(theta),sol(1));
% Plot results 
% Plot Pressure vs. Change in Theta
fig1 = figure('Name','Pressure vs. Change in Theta'); 
ax1 = axes(fig1); 
plot(ax1,theta,sol(:,1),'r-',MarkerSize=8);
hold on;
grid on; 
xlabel("Theta (radians)","FontSize",14);
ylabel("Pressure (Pa)","FontSize",14);
title("Pressure vs. Change in Theta","FontSize",14);
% Plot Volume vs. Change in Theta
fig2 = figure('Name','Volume vs. Change in Theta'); 
ax2 = axes(fig2); 
plot(ax2,theta,Vol1,'ro',MarkerSize=8);
hold on; 
grid on;
xlabel("Theta (radians)","FontSize",14);
ylabel("Volume (m^3)","FontSize",14);
title("Volume vs. Change in Theta","FontSize",14);
% Plot Work Out vs. Change in Theta
fig3 = figure('Name','Work Out vs. Change in Theta'); 
ax3 = axes(fig3); 
plot(ax3,theta,sol(:,2),'r-',MarkerSize=8);
hold on;
grid on;
xlabel("Theta (radians)","FontSize",14);
ylabel("Work Out (J)","FontSize",14);
title("Work Out vs. Change in Theta","FontSize",14);
% Plot Temperature vs. Change in Theta
fig4 = figure('Name','Temperature vs. Change in Theta'); 
ax4 = axes(fig4); 
plot(ax4,theta,Temp1(:,1),'ro',MarkerSize=8);
hold on; 
grid on; 
xlabel("Theta (radians)","FontSize",14);
ylabel("Temperature (K)","FontSize",14);
title("Temperature vs. Change in Theta","FontSize",14);

% Problem Number: 2 a)
% Update C and h values 
C = 0.8; % Empirical proportionality constant (s^-1)
h = 50; % Convective heat transfer parameter (W/m^2K)
% Solve for Pressure and Work In using ode45
[theta2,sol2] = ode45(@myeqs,theta_domain,y0);
% Solve for Volume (m^3)
[Vol2,dVol2] = myv(theta2);
% Solve for Temperature (K)
Temp2 = myt(Vol2,mym(theta2),sol2(1));
% Plot Results 
% Plot Pressure vs. Change in Theta
plot(ax1,theta2,sol2(:,1),'g-',MarkerSize=8);
hold on; 
% Plot Volume vs. Change in Theta
plot(ax2,theta2,Vol2,'gx',MarkerSize=8);
hold on;
% Plot Work Out vs. Change in Theta
plot(ax3,theta2,sol2(:,2),'g-',MarkerSize=8);
hold on; 
% Plot Temperature vs. Change in Theta
plot(ax4,theta2,Temp2(:,1),'gx',MarkerSize=8);
hold on; 

% Problem Number: 2 b)
% Update w value 
w = 100; % Engine speed (rad/s)
% Solve for Pressure and Work In using ode45
[theta3,sol3] = ode45(@myeqs,theta_domain,y0);
% Solve for Volume (m^3)
[Vol3,dVol3] = myv(theta3);
% Solve for Temperature (K)
Temp3 = myt(Vol3,mym(theta3),sol3(1));
% Plot Results 
% Plot Pressure vs. Change in Theta
plot(ax1,theta3,sol3(:,1),'b-',MarkerSize=8);
legend(ax1,'Pressure for Test 1','Pressure for Test 2 part a)','Pressure for Test 2 part b)','Location','best');
% Plot Volume vs. Change in Theta
plot(ax2,theta3,Vol3,'b*',MarkerSize=8);
legend(ax2,'Volume for Test 1','Volume for Test 2 part a)','Volume for Test 2 part b)','Location','best');
% Plot Work Out vs. Change in Theta
plot(ax3,theta3,sol3(:,2),'b-',MarkerSize=8);
legend(ax3,'Work Out for Test 1','Work Out for Test 2 part a)','Work Out for Test 2 part b)','Location','best');
% Plot Temperature vs. Change in Theta
plot(ax4,theta3,Temp3(:,1),'b*',MarkerSize=8);
legend(ax4,'Temperature for Test 1','Temperature for Test 2 part a)','Temperature for Test 2 part b)','Location','best');

% Functions Used %
    function derv = myeqs(theta,y)
        % y(1) = pressure
        [Vol,dVol] = myv(theta); % Values for Volume and the Derivative of Volume
        Mass = mym(theta); % Value for Mass
        dX = myx(theta); % Value for the derivative of the instantaneous fraction of total heat added to the gas
        Temp = myt(Vol,Mass,y(1)); % Value for Temperature 
        dP = myp(y(1),Vol,dVol,Mass,dX,Temp); % Value for Pressure 
        dW_out = y(1)*dVol; % Value for the Work Out 
        derv = [dP;dW_out];
    end
    function [V,dv_theta] = myv(theta) % Function that computes the volume and derivate of volume
        v = sym('v');
        a = sym('a');
        v(a) = V0*(1+((r-1)/(2*epsilon))*(1+(epsilon*(1-cos(a)))-(1-((epsilon)^2*(sin(a))^2))^(1/2))); % Volume Function
        V = double(v(theta)); % Value of Volume
        dv = diff(v,a); % Derivative of Volume
        dv_theta = double(dv(theta)); % Value of the Derivative of Volume 
    end
    function M = mym(theta) % Function that computes the mass
        m = sym('m');
        c = sym('c');
        m(c) = M0*exp(-(C/w)*(c-pi)); % Function for mass
        M = double(m(theta)); % Value for mass 
    end
    function dX_theta = myx(theta) % Function for computing the derivative of the instantaneous fraction of total heat added to the gas
        if (theta >= pi) && (theta < theta_s)
            X = sym('X');
            d = sym('d');
            X(d) = 0; % Function for the instantaneous fraction of total heat added to the gas
            dX = diff(X,d); % Derivative of the the instantaneous fraction of total heat added to the gas
            dX_theta = double(dX(theta)); % Value for the derivative of the instantaneous fraction of total heat added to the gas
        elseif (theta >= theta_s) && (theta <= (theta_s + theta_b))
            X = sym('X');
            d = sym('d');
            X(d) = (1/2)*(1-cos((pi*(d-theta_s))/theta_b)); % Function for the instantaneous fraction of total heat added to the gas
            dX = diff(X,d); % Derivative of the the instantaneous fraction of total heat added to the gas
            dX_theta = double(dX(theta)); % Value for the derivative of the instantaneous fraction of total heat added to the gas
        else
            X = sym('X');
            d = sym('d');
            X(d) = 1; % Function for the instantaneous fraction of total heat added to the gas
            dX = diff(X,d); % Derivative of the the instantaneous fraction of total heat added to the gas
            dX_theta = double(dX(theta)); % Value for the derivative of the instantaneous fraction of total heat added to the gas
        end
    end
    function t = myt(V,M,P) % Function for computing temperature 
        t = (P*V)/(M*R); % Value for temperature (K) 
    end
    function m_dot = myL(M) % Function for computing the instantaneous mass loss rate
        m_dot = C*M; % Value for the instantaneous mass loss rate (kg/s)
    end
    function A = myA(V) % Function for computing the instantaneous surface area 
        A = (4*V)/b; % Value for the instantaneous surface area (m^2)
    end
    function dp = myp(P,V,dV,M,dx,T) % Function for computing the derivative of pressure 
        dp = (-gamma*(P/V)*dV)+((gamma-1)*((M*q_in)/V)*dx)-(((gamma-1)*h*myA(V)*(T-Tw))/(V*w))-(gamma*(myL(M)/(M*w))*P); % Value for derivative of pressure
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%