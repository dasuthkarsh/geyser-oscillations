format long;
clear;
%close all;

addpath ./XSteam_Matlab_v2.6/;
addpath ./Steam/
addpath ..;

par.sb = 80; % cross-section of bubble trap in m^2
par.sc = 1; % cross-section of geyser column in m^2
par.H = 7; % height of bubble trap in m
par.L = 0; % length of lateral connector in m
par.sl = 1; %cross-section of lateral connector in m^2
par.xbar = par.H-0.004539970908913; %mean water height in bubble trap in m, must be in (0,H)
par.ybar = 15;
par.delxy=par.ybar-par.xbar;

par.g = 9.81;% gravitational acceleration in m*s^-2
par.rho = 1000; % water density in kg*m^-3
par.gamma=7/5; % adiabatic exponent for diatomic ideal gas - unitless
par.alpha = 5/2; %diatomic ideal gas constant - unitless
par.Pa0 =1e5; %atmospheric pressure at equilibrium in Pa

numPar.v0= 0; %initial velocity in m/s
numPar.tf = 5;
numPar.x0 = -1/80;

%% Make a lookup table with pre-computed thermodynamic properties.
par.Vol_0 = par.sb*(par.H-par.xbar); %initial (equilibrium) volume in m^3
P_0 = par.rho*par.g*(par.delxy)+par.Pa0; %initial pressure in (Pascal)
par.m=par.Vol_0 / XSteam('vV_p', P_0/1e5); %vapor mass in kg

xv = 1 - 1e-2;
sv=XSteam('sV_p',P_0/1e5)*1e3; %specific entropy J*kg^-1*degC^-1
sl=XSteam('sL_p',P_0/1e5)*1e3; %specific entropy J*kg^-1*degC^-1
s = xv*sv + (1-xv)*sl; % gives specific entropy of water+vapor mixture in J/kg/C
V0 = XSteam('V_ps',P_0/1e5,s/1e3); % volume per unit mass
par.m = par.Vol_0 / V0;

Pmax=[9e5,2e6];
Pmin=[9e4,1.6e5];

for j=1:length(Pmin)
    P=linspace((P_0-Pmin(j))/1e5,(P_0+Pmax(j))/1e5,101); %Range of pressure for lookup table (bar)
    vV=zeros(1,length(P));  % vapor volume (m^3)
    du_dp = zeros(1,length(P)); %dU/dP in (kJ/K/Pa)
    dv_dp = zeros(1,length(P)); %dV/dP in (m^3/Pa)
    xS = zeros(1,length(P));
    T = zeros(1,length(P));
    for i=1:length(P)
        vV(i) = XSteam('v_ps',P(i),s/1e3)*par.m; %volume in m^3
        xS(i) = XSteam('x_ps',P(i),s/1e3); %vapor fraction (unitless)
        T(i) = XSteam('Tsat_p',P(i)); % temperature in degree C
        delta = 1e-6;%bar
        du_dp(i) = (XSteam('u_ps',P(i)+delta,s/1e3)-XSteam('u_ps',P(i)-delta,s/1e3))*1e3*par.m/(2*delta*1e5);% enthalpy in mks units, pressure in Pa
        dv_dp(i) = (XSteam('v_ps',P(i)+delta,s/1e3)-XSteam('v_ps',P(i)-delta,s/1e3))*par.m/(2*delta*1e5);% m^3/Pa
    end
    [~,i] = sort(vV);
    par.Fdudp = griddedInterpolant(vV(i),du_dp(i),'linear','none');
    par.FdPdv = griddedInterpolant(vV(i),1./dv_dp(i),'linear','none');
    if (j==1)
        fprintf("Steam tables computed for pressures ranging from %f bar to %f bar \n", (P_0-Pmin(j))/1e5,(P_0+Pmax(j))/1e5);
        [freq,~]= frequency_calculator(par);
        fprintf("steam_model would predict a frequency of %f \n", freq);
    elseif (j==2)
        fprintf("Steam tables computed for pressures ranging from %f bar to %f bar\n", (P_0-Pmin(j))/1e5,(P_0+Pmax(j))/1e5);
        [freq,~]= frequency_calculator(par);
        fprintf("application_to_ofg would predict a frequency of %f \n", freq);
    end
end