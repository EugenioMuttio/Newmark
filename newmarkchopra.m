
clc
clear all
format long
%Time plot

timeplot=10;

%Newmark Method
%Initial Conditions
u0=0;
v0=1;
a0=0;

%Parameter
gamma=0.5;
beta=0.25;
w=0.5;
T=1/w;
%time
tf=20*T;
delta_t=0.1;


%Exact Solution
sol_ex=@(t) v0*sin(w*t)/w;

u_n1=[];
v_n1=[];
a_n1=[];
p_hat=[];

u_n1(1)=u0;
v_n1(1)=v0;
a_n1(1)=-w^2*u_n1(1);

it=delta_t;
n=tf/delta_t;
vt = linspace(0,tf,n+1);

k1=1/(beta*delta_t^2);
k2=1/(beta*delta_t);
k3=(1/(2*beta)-1);
k_hat=w^2+k1;
i=2; %index 

while it<=tf+delta_t
    p_hat(i)=k1*u_n1(i-1)+k2*v_n1(i-1)+k3*a_n1(i-1);
    
    u_n1(i)=p_hat(i)/k_hat;
    v_n1(i)=gamma/(beta*delta_t)*(u_n1(i)-u_n1(i-1))+(1-gamma/beta)*v_n1(i-1)+delta_t*(1-gamma/(2*beta))*a_n1(i-1);
    a_n1(i)=1/(delta_t^2*beta)*(u_n1(i)-u_n1(i-1))-1/(beta*delta_t)*v_n1(i-1)-(1/(2*beta)-1)*a_n1(i-1);
    it=it+delta_t;
    i=i+1;
end



%Newmark Displacement Plot
figure(1)
hold on
plot(vt,u_n1, '-r');
plot(vt,v_n1, '-k');
tit='Solution Comparison beta='+string(beta)+' $\Delta t=$'+string(delta_t);
title(tit,'Interpreter','latex','FontSize',17);
xlabel('$t$','Interpreter','latex','FontSize',17);
ylabel('$u$','Interpreter','latex','FontSize',17);

%Exact Displacement Plot
n_ex=tf/0.001;
vt_ex = linspace(0,tf,n+1);
plot(vt_ex,sol_ex(vt_ex),'-b');
grid on
legend('Newmark','Velocity','Exact')
hold off

