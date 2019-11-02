clc
clear all
format long
%Time plot

timeplot=10;

%Newmark Method
%Initial Conditions
u0=0;
v0=5;
a0=0;

%Parameter
gamma=0.5;
beta=[0,1/12,1/6,1/4,1/3];
w=0.5;
T=1/w;
%time
tf=20*T;


for ib=1:length(beta)
    delta_t=[0.001,0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.55,0.60,0.7,0.8,0.9,1];

    deltaT=[];
    KE=[];

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
        
        maxvel=[];

        u_n1(1)=u0;
        v_n1(1)=v0;
        a_n1(1)=a0;
             

        it=delta_t(dti);
        n=tf/delta_t(dti);
        vt = linspace(0,tf,n+1);

        i=2; %index 
        iP=1; % Newmark Period Counter
        iP2=0; % Exact Period Counter

        while it<=tf+delta_t

            u_n1(i)=u_n1(i-1) + delta_t(dti)*v_n1(i-1)+0.5*delta_t(dti)^2*(1-2*beta(ib))*a_n1(i-1);
            v_n1(i)=v_n1(i-1)+delta_t(dti)*(1-gamma)*a_n1(i-1);
            a_n1(i)=-w^2*u_n1(i)/(1+w^2*beta(ib)*delta_t(dti)^2);

            %Newmark Period
            if sign(u_n1(i))~=sign(u_n1(i-1)) 
                t1=it-delta_t(dti);
                t2=it;
                vP(iP)=((t2-t1)/(u_n1(i)-u_n1(i-1)))*(0-u_n1(i-1))+t1;

                %Exact 
                P_ex(iP)=pi*iP2/w;

                iP=iP+1;
                iP2=iP2+1;

            end

 
            it=it+delta_t(dti);
            i=i+1;

        end

        %Error Periodo en 3 indice
        deltaT(dti)=(vP(3)-P_ex(3));
        deltaTT(dti)=deltaT(dti)/T;
        delta_tT(dti)=delta_t(dti)/T;
        
        %Kinetic Energy
        maxvel(dti)=max(abs(v_n1(:)));
        KE(dti)=0.5*maxvel(dti)^2;
        
        if dti==timeplot
            %Newmark Displacement Plot
            figure(ib)
            hold on
            plot(vt,u_n1, '-r');
            plot(vt,v_n1, '-k');
            tit='Solution Comparison beta='+string(beta(ib))+' $\Delta t=$'+string(delta_t(timeplot));
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
        end

    end
    
    
    
    %Error Plots
    figure(6)
    hold on
    plot(delta_tT,deltaTT);
    title('Error Plot','Interpreter','latex','FontSize',17);
    legend('beta=0','beta=1/12','beta=1/6','beta=1/4','beta=1/3');
    xlabel('$\Delta t / T$','Interpreter','latex','FontSize',17);
    ylabel('$\Delta T / T$','Interpreter','latex','FontSize',17);
    grid on
    hold off 
    
    %Kinetic Plots
    figure(7)
    hold on
    plot(log10(delta_t),log10(KE));
    title('Kinetic Energy Plot','Interpreter','latex','FontSize',17);
    legend('beta=0','beta=1/12','beta=1/6','beta=1/4','beta=1/3');
    xlabel('$log_{10}(\Delta t)$','Interpreter','latex','FontSize',17);
    ylabel('$log_{10}(K)$','Interpreter','latex','FontSize',17);
    grid on
    hold off 
    

    
end








