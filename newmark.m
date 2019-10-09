%Newmark Method

%Initial Conditions
 

u0=0;
v0=2;
a0=0;

%Parameter
gamma=0.5;
beta=1/12;

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

u_n1(1)=u0;
v_n1(1)=v0;
a_n1(1)=a0;

u_ex=[];
vt=[];
u_n1_b=[];
deltaT=[];
deltaT_b=[];

%Exact Solution

sol_ex=@(t) v0*sin(w*t)/w;

n=tf/delta_t;
vt = linspace(0,tf,n+1);
nT=tf/T;
vP = [];
i=2;
iP=1;
iter=0;

%for b=1:length(delta_t)
    while iter<=tf   

        u_n1(i)=un + delta_t*vn+0.5*delta_t^2*(1-2*beta)*an;
        v_n1(i)=vn+delta_t*(1-gamma)*an;
        a_n1(i)=-w^2*u_n1(i)/(1+w^2*beta*delta_t^2);
        
        if sign(u_n1(i))~=sign(un) 
            t1=iter;
            t2=iter+delta_t;
            vP(iP)=((t2-t1)/(u_n1(i)-un))*(0-un)+t1;
            iP=iP+1;
        end

        un=u_n1(i);
        vn=v_n1(i);
        an=a_n1(i);
        

        iter=iter+delta_t;
        i=i+1;
        
    end

%end



plot(vt,u_n1, '-r');
hold on
%plot(vt,sol_ex(vt),'-b');
plot(vt,v_n1, '-k');
plot(vP,sol_ex(0), 'og');
grid on
hold off



legend('Newmark','Exact','vel')

