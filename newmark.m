clc
clear all 
%Newmark Method

%Initial Conditions

u0=0;
v0=2;
a0=0;

%Parameter
gamma=0.5;
beta=1/12;

w=1;
T=1/w;

%time
tf=100*T;
delta_t=[0.05,0.1,0.2,0.5];

deltaT=[];

u_n1(1)=u0;
v_n1(1)=v0;
a_n1(1)=a0;

u_ex=[];
vt=[];


%Exact Solution
sol_ex=@(t) v0*sin(w*t)/w;







for dti=1:length(delta_t)
    
    u_n1=[];
    v_n1=[];
    a_n1=[];
    vP = [];
    P_ex=[];
    P_ex(1)=0;

    
    un=u0;
    vn=v0;
    an=a0;
    
    it=delta_t(dti);
    n=tf/delta_t(dti);
    vt = linspace(0,tf,n+1);
    
    i=2;
    iP=1;
    iP2=0;
    
    while it<=tf

        u_n1(i)=un + delta_t(dti)*vn+0.5*delta_t(dti)^2*(1-2*beta)*an;
        v_n1(i)=vn+delta_t(dti)*(1-gamma)*an;
        a_n1(i)=-w^2*u_n1(i)/(1+w^2*beta*delta_t(dti)^2);
        
        %Newmark Period
        if sign(u_n1(i))~=sign(un) 
            t1=it-delta_t(dti);
            t2=it;
            vP(iP)=((t2-t1)/(u_n1(i)-un))*(0-un)+t1;
            
            %Exact 
            P_ex(iP)=pi*iP2/w;
            
            iP=iP+1;
            iP2=iP2+1;
            
        end

        un=u_n1(i);
        vn=v_n1(i);
        an=a_n1(i);
        

        it=it+delta_t(dti);
        i=i+1;
        
    end
    
    %Error Periodo en 3 indice
    
    deltaT(dti)=(vP(3)-P_ex(3));
    
    
    
end


%Newmark Plot
plot(vt,u_n1, '-r');
hold on
%Exact Solution
plot(vt,sol_ex(vt),'-b');
%Velocity Plot
%plot(vt,v_n1, '-k');

plot(vP,sol_ex(0), 'og');
%plot(P_ex,sol_ex(P_ex), '+k');
grid on
legend('Newmark','Exact','vel')
hold off





