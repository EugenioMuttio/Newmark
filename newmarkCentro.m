clc
clear all
close all
format long


%Newmark Method

%Earthquake data
T=readtable('elcentro.dat');
delta_t=T{3,1}-T{2,1};
accsize=size(T);
accsize=accsize(1);

grav=980.665;

%Initial Conditions
u0=0;
v0=0;
a0=0;

u_max=[];
v_max=[];
a_max=[];

%Newmark Parameters
gamma=0.5;
beta=0;

%Dynamic Sistem Parameter
xi=0.1;
NTs=10;
NT=linspace(0.1,10,NTs); %Natural Period
w=[];

%Plots by chosen frequency
wp=5;

for i=1:NTs
    w(i)=1/NT(i);
end

for wi=1:NTs
    u_n1=[];
    v_n1=[];
    a_n1=[];

    u_n1(1)=u0;
    v_n1(1)=v0;
    a_n1(1)=a0;

    tf=T{accsize(1),1};
    n=tf/delta_t;
    vt = linspace(0,tf,n+1);

    i=2; %index 
    iP=1; % Newmark Period Counter
    it=delta_t;
    while i<=accsize

        u_n1(i)=u_n1(i-1) + delta_t*v_n1(i-1)+0.5*delta_t^2*(1-2*beta)*a_n1(i-1);
        v_n1(i)=v_n1(i-1)+delta_t*(1-gamma)*a_n1(i-1);
        a_n1(i)=(T{i,2}*grav-2*xi*w(wi)*v_n1(i)-w(wi)^2*u_n1(i))/(1+gamma*2*xi*w(wi)*delta_t+w(wi)^2*beta*delta_t^2);

        it=it+delta_t;
        i=i+1;
    end
    
    %Maximum Values
    u_max(wi)=max(abs(u_n1));
    v_max(wi)=max(abs(v_n1));
    a_max(wi)=max(abs(a_n1));
    
    if wi==wp
        %Newmark Displacement Plot
        figure(1)
        plot(vt,u_n1, '-r');
        title('Displacement Plot','Interpreter','latex','FontSize',19);
        xlabel('$t$','Interpreter','latex','FontSize',19);
        ylabel('$u(t)$','Interpreter','latex','FontSize',19);
        grid on

        %Newmark Velocity Plot
        figure(2)
        plot(vt,v_n1, '-b');
        title('Velocity Plot','Interpreter','latex','FontSize',19);
        xlabel('$t$','Interpreter','latex','FontSize',19);
        ylabel('$v(t)$','Interpreter','latex','FontSize',19);
        grid on

        %Newmark Acceleration Plot
        figure(3)
        plot(vt,a_n1, '-k');
        title('Acceleration Plot','Interpreter','latex','FontSize',19);
        xlabel('$t$','Interpreter','latex','FontSize',19);
        ylabel('$a(t)$','Interpreter','latex','FontSize',19);
        grid on

    
    end
    
    
end


% Accelerogram Plot
figure(4)
plot(T{:,1},T{:,2}*grav, '-b');
title('El centro Plot','Interpreter','latex','FontSize',19);
xlabel('$t [s]$','Interpreter','latex','FontSize',19);
ylabel('$a_{g}(t) [cm/s^2]$','Interpreter','latex','FontSize',19);
grid on

%Plot uMax-Periods
figure(5)
plot(NT,u_max, '-*r');
title('Max Displacement - Natural Period','Interpreter','latex','FontSize',19);
xlabel('$t$','Interpreter','latex','FontSize',19);
ylabel('$u_{max}$','Interpreter','latex','FontSize',19);
grid on

%Plot vMax-Periods
figure(6)
plot(NT,v_max, '-*r');
title('Max Velocity - Natural Period','Interpreter','latex','FontSize',19);
xlabel('$t$','Interpreter','latex','FontSize',19);
ylabel('$v_{max}$','Interpreter','latex','FontSize',19);
grid on

%Plot uMax-Periods
figure(7)
plot(NT,a_max, '-*r');
title('Max Acceleration - Natural Period','Interpreter','latex','FontSize',19);
xlabel('$t$','Interpreter','latex','FontSize',19);
ylabel('$a_{max}$','Interpreter','latex','FontSize',19);
grid on




    









