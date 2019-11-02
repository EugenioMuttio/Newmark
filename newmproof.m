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
delta_t=0.001;


%Exact Solution
sol_ex=@(t) v0*sin(w*t)/w;

u_n1=[];
v_n1=[];
a_n1=[];

u_n1(1)=u0;
v_n1(1)=v0;
a_n1(1)=a0;

it=delta_t;
n=tf/delta_t;
vt = linspace(0,tf,n+1);

i=2; %index 

while it<=tf

    u_n1(i)=u_n1(i-1) + delta_t*v_n1(i-1)+0.5*delta_t^2*(1-2*beta)*a_n1(i-1);
    v_n1(i)=v_n1(i-1)+delta_t*(1-gamma)*a_n1(i-1);
    a_n1(i)=-w^2*u_n1(i)/(1+w^2*beta*delta_t^2);
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







