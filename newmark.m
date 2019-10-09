%Newmark Method

%Initial Conditions
 

u0=0;
v0=2;
a0=0;

%Parameter
gamma=0.5;
beta=0;

w=5;
T=1/w;

%time
ti=1;
tf=100*T;
delta_t=0.1;

un=u0;
vn=v0;
an=a0;

u_n1=[];
v_n1=[];
a_n1=[];
u_ex=[];
vt=[];

%Exact Solution

sol_ex=@(t) v0*sin(w*t)/w;

for i=ti:tf    
    
    u_n1(i)=un + delta_t*vn+0.5*delta_t^2*(1-2*beta)*an;
    v_n1(i)=vn+delta_t*(1-gamma)*an;
    a_n1(i)=-w^2*u_n1(i)/(1+w^2*beta*delta_t^2);
    
    un=u_n1(i);
    vn=v_n1(i);
    an=a_n1(i);
    
end

n=tf/delta_t;
vt = linspace(0,tf,n);

plot(u_n1);
hold on
plot(vt,sol_ex(vt));
plot(vt,2/5*sin(5*vt));

