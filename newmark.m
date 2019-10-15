clc
clear all
%Newmark Method
%Initial Conditions
u0=0;
v0=2;
a0=0;

%Parameter
gamma=0.5;
beta=[0,1/12,1/6,1/4,1/3];
w=1;
T=1/w;
%time
tf=100*T;

for ib=1:length(beta)
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

        i=2; %index 
        iP=1; % Newmark Period Counter
        iP2=0; % Exact Period Counter

        while it<=tf

            u_n1(i)=un + delta_t(dti)*vn+0.5*delta_t(dti)^2*(1-2*beta(ib))*an;
            v_n1(i)=vn+delta_t(dti)*(1-gamma)*an;
            a_n1(i)=-w^2*u_n1(i)/(1+w^2*beta(ib)*delta_t(dti)^2);

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
        deltaTT(dti)=deltaT(dti);
        delta_tT(dti)=delta_t(dti);

    end

    %Error Plots
    figure(2)
    hold on
    plot(delta_tT,deltaTT);
    legend('beta=0','beta=1/12','beta=1/6','beta=1/4','beta=1/3');
    grid on
    hold off 
end


%Newmark Plot
figure(1)
plot(vt,u_n1, '-r');
hold on
%Exact Solution
n_ex=tf/0.001;
vt_ex = linspace(0,tf,n+1);
plot(vt_ex,sol_ex(vt_ex),'-b');
grid on
legend('Newmark','Exact')
hold off



