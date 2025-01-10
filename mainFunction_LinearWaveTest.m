%%%%%%%%%%%%%%%%%%%%%%%copy right @ Longhuan Zhu 2020%%%%%%%%%%%%%%%%%%%%%%
%---------------------------version 2020.01--------------------------------
%
% This is a MATLAB transcript to calculate the dynamics of flexible blade in
% unsteady flow (wave with/without current). This code is developed by 
% Dr. Longhuan Zhu to analyze the asymmetric motion of vegetation in waves 
% (Zhu et al., 2020). The code was validated by laboratory experiments in
% Luhar and Nepf (2016) and Zeller et al. (2014).
%  
% Recommended Citation:
% Zhu, L., Zou, Q., Huguenard, K., & Fredriksson, D. W. (2020).
% Mechanisms for the Asymmetric Motion of Submerged Aquatic Vegetation in
% Waves: A Consistent-Mass Cable Model. JGR: Oceans, 125(2).
% https://doi.org/10.1029/2019JC015517
%
% Contact longhuan.zhu@unh.edu in case of any questions.
%
% reference:
% [1] Luhar, M., &Nepf, H. M. (2016). Wave-induced dynamics of flexible blades.
% Journal ofFluids and Structures, 61, 20–41.
% https://doi.org/10. 1016/j.jfluidstructs.2015.11.007
%
% [2] Zeller, R. B., Weitzman, J. S., Abbett, M. E., Zarama, F. J., 
% Fringer, O. B., &Koseff, J. R. (2014). Improved parameterization of
% seagrass blade dynamics and wave attenuation based on numerical and 
% laboratory experiments. Limnology and Oceanography, 59(1), 251–266. 
% https:// doi.org/10.4319/lo.2014.59.1.0251
%
%%%%%%%%%%%%%%%%%%%%%%%copy right @ Longhuan Zhu 2020%%%%%%%%%%%%%%%%%%%%%%
%
%% coordinate system
% Fixed global Cartesian reference frame (x, z):
% Origin at the blade base, x in the horizontal direction, z in vertical directions
% The flow field is described by the horizontal and vertical components, 
% U(x, z, t) and W(x, z, t) with the fixed global coorindate system
% 
% Local Lagrangian coordinate system (t, n):
% On local blade segment, t representing the blade-tangential direction and
% n as the blade-normal direction. The variables for the blade dynamics is 
% defined with the local Lagrangian coorindate system [T,Q,u,v,Omega,phi]
% with T is the effective tension in the blade-tangential direction, Q the
% shear force in the blade-normal direction, u the velocity component in
% the blade-tangential direction, v the velocity component in the
% blade-normal direction, Omega the blade curvature, phi the bending angle.
% Note that the angle here is defiened as the angle from x to t, while the
% angle in Zhu et al. (2020) is defined as the angle from vertical to t for
% the convenience of asymmetric motion analysis. [T,Q,u,v,Omega,phi] are
% the varaibles being solved.

%% input variables in metric units (SI) and examples
% %constant parameters
% g=9.81;             % [m s^-2] gravatational acceleration
% tol=1e-8;           % Tolerance    
%
% %Structural Parameters (geometrical and mechanical properties of blade)
% length=0.2;         % [m] structure length(vegetation height) 20cm
% width=0.02;         % [m] structure width: this case for rectanglur cross section 2cm   
% thick=0.0019;       % [m] structure thickness
% rhou_structure=670; % [kg m^-3] structure mass density
% E=5e5;              % [Pa] structure elastic modulus (Young's modulus)
% 
% %Flow Parameters (waves and currents)
% rhou_water=1000;    % [kg m^-3] water mass density
% Depth_water=0.3;    % [m] water depth
% Height_wave=0.08;   % [m] wave height
% T_wave=2;           % [m] wave period
% Uc=0;               % [m/s] current speed
% 
% %hydrodynamic coefficients
% Ca=1;               % added mass coeff.
% Cdn=nan;%1.95;      % drag coeff. (if nan, Cdn and Ca will be calculated 
%                     %             with the formulas 12 and 13 in Zhu et al., JGR2020)
% Cdt=nan;%0.01;      % friction coeff. (if nan, Cdt will be calcualted 
%                     %             with the formula 14 in Zhu et al., JGR2020)
% 
% %control options
% n=20;               % number of segments   
% Dt=0.01;            % time step
% num_T_period=4;     % number of periods to simulate (total simulation time)
%                     %  It should be large enough to get into steady state
% wave_option=0;      %0: measured water partical velocity
%                     %1: linear wave theory
%                     %2: user defined
% 

%% output
% (DisX, DisZ): [m] blade position coorindates
% [T,Q,u,v,Omega,phi]
% T: [N] effective tension in the blade-tangential direction
% Q: [N] shear force in the blade-normal direction
% u: [m/s] velocity component in the blade-tangential direction
% v: [m/s] the velocity component in the blade-normal direction
% Omega: [m^-1] blade curvature
% phi: [radians] blade bending angle.
%
% optional output
% DragT: [N/m] drag (friction) in the blade-tangential direction
% DragN: [N/m] drag (normal drag) in the blade-normal direction
% AMN: [N/m] added mass force 
% (VBT, VBN): [N/m] virtual buoyance (Froude–Krylov force)

%% input
%constant parameters
g=9.81;             % [m s^-2] gravatational acceleration
tol=1e-8;           % Tolerance    

%Structural Parameters (geometrical and mechanical properties of blade)
length=0.2;         % [m] structure length(vegetation height) 20cm
width=0.02;         % [m] structure width: this case for rectanglur cross section 2cm   
thick=0.0019;       % [m] structure thickness
rhou_structure=670; % [kg m^-3] structure mass density
E=5e5;              % [Pa] structure elastic modulus (Young's modulus)

%Flow Parameters (waves and currents)
rhou_water=1000;    % [kg m^-3] water mass density
Depth_water=0.3;    % [m] water depth
Height_wave=0.08;   % [m] wave height
T_wave=2;           % [m] wave period
Uc=0;               % [m/s] current speed

%hydrodynamic coefficients
Ca=1;               % added mass coeff.
Cdn=nan;%1.95;      % drag coeff. (if nan, Cdn and Ca will be calculated 
                    %             with the formulas 12 and 13 in Zhu et al., JGR2020)
Cdt=nan;%0.01;      % friction coeff. (if nan, Cdt will be calcualted 
                    %             with the formula 14 in Zhu et al., JGR2020)

%control options
n=20;               % number of segments   
Dt=0.01;            % time step
num_T_period=4;     % number of periods to simulate (total simulation time)
                    %  It should be large enough to get into steady state
wave_option=1;      %0: measured water partical velocity from Luhar
                    %1: linear wave theory
                    %2: user defined in line 138-143
    
%% %%%%%%%%%%%%%% start to calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% initial calculation
%structure
Area_cross=width*thick;             %cross section area
I=width*thick^3/12;                 %second momentum for rectangular cross
m_unit=rhou_structure*Area_cross;   %mass for unit length
%         ma_unit=Ca*rhou_water*pi/4*width^2; %added mass for unit length
ma_unit=Ca*rhou_water*pi/4*width^2; %width*thick;%added mass for unit length
w0_unit=(rhou_structure-rhou_water)*Area_cross*g;                   %wet weight for unit length

%waves        
Omega_wave=2*pi/T_wave;             %angular frequency        
k_wave_number=wave_num(Omega_wave,Depth_water);%wave number (linear wave theory)

%linear wave function
if wave_option < 1
elseif wave_option < 2
    Ux_wave=@(x,z,t) Height_wave/2*Omega_wave*cosh(k_wave_number.*z)/sinh(k_wave_number*Depth_water).*cos(k_wave_number.*x-Omega_wave.*t) + Uc;
    Uz_wave=@(x,z,t) Height_wave/2*Omega_wave*sinh(k_wave_number.*z)/sinh(k_wave_number*Depth_water).*sin(k_wave_number.*x-Omega_wave.*t) + Uc;
    Ux_max=@(x,z,t) Height_wave/2*Omega_wave*cosh(k_wave_number.*z)/sinh(k_wave_number*Depth_water) + abs(Uc);
elseif wave_option < 3
    % user define
    % Ux_wave=@(x,z,t) Ux
    % Uz_wave=@(x,z,t) Uz
    % Ux_max=@(x,z,t) Ux max
end

%control
T_total=num_T_period*T_wave;    %total time
num_time=T_total/Dt;
t=0:Dt:T_total;                 %time

%% Main Procedure
%initialize
Ds=length/n;                    %segment length
m(:,1)=zeros(n,1)+m_unit;       %segment mass
ma_seg=ma_unit;                   %segment added mass
w0(:,1)=zeros(n,1)+w0_unit;%*Ds;  %segment weight

% unknow vector: X=[T,Q,u,v,Omega,phi]
T=zeros(n+1,num_time+1);        %tension
Q=zeros(n+1,num_time+1);        %shear
u=zeros(n+1,num_time+1);        %tangent velocity
v=zeros(n+1,num_time+1);        %normal velocity
Omega=zeros(n+1,num_time+1);    %curvature
phi=zeros(n+1,num_time+1);      %angular from horizontal

%time=0;
T(n+1,1)=0;                     %-rhou_water*g*(Depth_water-length)*Area_cross;            %0;
for i=n+1:-1:2
    T(i-1,1)=T(i,1)+w0(i-1)*Ds;%+rhou_water*g*(Depth_water-(i-2)*Ds)*Area_cross; %pretension
end
phi(:,1)=zeros(n+1,1)+pi/2;     %vertical pi/2

   
%solution definition X=[T,Q,u,v,Omega,phi]
X=zeros(6*(n+1),1);
for k=1:n+1                     %kth segment
    X(6*(k-1)+1)=T(k,1);
    X(6*(k-1)+2)=Q(k,1);
    X(6*(k-1)+3)=u(k,1);
    X(6*(k-1)+4)=v(k,1);
    X(6*(k-1)+5)=Omega(k,1);
    X(6*(k-1)+6)=phi(k,1);
end

%postion    
DisX=zeros(n+1,num_time+1);     %x displacement for each point
DisZ=zeros(n+1,num_time+1);     %z displacement for each point
DisZ(1,1)=0;
for k=2:n+1
        DisZ(k,1)=DisZ(k-1,1)+Ds;
end
    
%structure length--check
Length_curvature=zeros(num_time+1,1);
Length_curvature(1)=length;

u_flow=zeros(n+1,num_time+1);       %water partical velocities: tangent direction
v_flow=zeros(n+1,num_time+1);       %normal direction
%     global Umax
Umax=zeros(n+1,num_time+1);

%select wave theory
if wave_option<1
    [u_flow(:,1),v_flow(:,1),Umax(:,1)]=wave_nearVege(DisX(:,1),DisZ(:,1),t(1),k_wave_number,Omega_wave);
else
    u_flow(:,1)=Ux_wave(DisX(:,1),DisZ(:,1),t(1));
    v_flow(:,1)=Uz_wave(DisX(:,1),DisZ(:,1),t(1));
    Umax(:,1)=Ux_max(DisX(:,1),DisZ(:,1),t(1)); 
end

%% time after 0
for j=2:num_time+1

    X0=X;

    %New position
    DisX(:,j)=DisX(:,j-1);
    DisZ(:,j)=DisZ(:,j-1);
    u_flow(:,j)=u_flow(:,j-1);
    v_flow(:,j)=v_flow(:,j-1);
    Umax(:,j)=Umax(:,j-1);

    %Substitude parameters into Governing Equations
    Eq_Govern_Done=@(X0,X,Ux_wave_old,Ux_wave_new,Uz_wave_old,Uz_wave_new,m,w0,Umax) Eq_Govern(X0,X,Ux_wave_old,Ux_wave_new,Uz_wave_old,Uz_wave_new,m,ma_seg,Dt,Ds,w0,rhou_water,Cdt,Cdn,width,Area_cross,E,I,thick,T_wave,Umax);

    Y=Y_value(Eq_Govern_Done,X0,X,u_flow(:,j-1),u_flow(:,j),v_flow(:,j-1),v_flow(:,j),n,m,w0,DisZ(n+1,j),Umax(:,j));
    % [Y,KC(:,j),~,DragT(:,j),DragN(:,j),MaT(:,j),MaN(:,j),BT(:,j),BN(:,j),TT(:,j),TN(:,j),VBT(:,j),WT(:,j),AMN(:,j),VBN(:,j),WN(:,j)]=Y_value(Eq_Govern_Done,X0,X,u_flow(:,j-1),u_flow(:,j),v_flow(:,j-1),v_flow(:,j),n,m,w0,DisZ(n+1,j),Umax(:,j));

    norm_Y=norm(Y);

    while norm_Y>tol
            
        %Form Jacobi
        J=J_value(@Jacob,Eq_Govern_Done,X0,X,u_flow(:,j-1),u_flow(:,j),v_flow(:,j-1),v_flow(:,j),n,m,w0,Umax(:,j));
        
        X=X-J\Y;

        %New position
        [DisX(:,j),DisZ(:,j),Length_curvature(j)]=Position(DisX(:,j-1),DisZ(:,j-1),X0,X,Dt,Ds,n);

        %select wave theory
        if wave_option<1
            [u_flow(:,j),v_flow(:,j),Umax(:,j)]=wave_nearVege(DisX(:,j),DisZ(:,j),t(j),k_wave_number,Omega_wave);
        else
            u_flow(:,j)=Ux_wave(DisX(:,j),DisZ(:,j),t(j));
            v_flow(:,j)=Uz_wave(DisX(:,j),DisZ(:,j),t(j));     
            Umax(:,j)=Ux_max(DisX(:,j),DisZ(:,j),t(j));         
        end
           
        Y=Y_value(Eq_Govern_Done,X0,X,u_flow(:,j-1),u_flow(:,j),v_flow(:,j-1),v_flow(:,j),n,m,w0,DisZ(n+1,j),Umax(:,j));
        % [Y,KC(:,j),~,DragT(:,j),DragN(:,j),MaT(:,j),MaN(:,j),BT(:,j),BN(:,j),TT(:,j),TN(:,j),VBT(:,j),WT(:,j),AMN(:,j),VBN(:,j),WN(:,j)]=Y_value(Eq_Govern_Done,X0,X,u_flow(:,j-1),u_flow(:,j),v_flow(:,j-1),v_flow(:,j),n,m,w0,DisZ(n+1,j),Umax(:,j));
           
        norm_Y=norm(Y);
    end

    %output results for jth time
    for k=1:n+1
        T(k,j)=X(6*(k-1)+1);
        Q(k,j)=X(6*(k-1)+2);
        u(k,j)=X(6*(k-1)+3);
        v(k,j)=X(6*(k-1)+4);
        Omega(k,j)=X(6*(k-1)+5);
        phi(k,j)=X(6*(k-1)+6);
    end
    
    %process report    
    fprintf('T_wave= %4.1fs; Time is %8.4fs\n',T_wave,t(j)) 
end

%% save data
file='test';
direction=['.\LinearWaveTestResults\',file];
mkdir(direction);
Filename=[file,'.mat'];
save([direction,'\',Filename]);
    
%Output 
%structural length check
AveSum=mean(Length_curvature);
AveSum2=sqrt(norm(Length_curvature)^2/(num_time+1));

%% plot the blade posture in the last wave period
figure(1)
clf;
T_num=T_wave/Dt;
DT=round(T_num/64);
%experiment data form Luhar and Nepf(2015)        
% Ex=load('./LuharData/exFoam.txt');
% Num=load('./LuharData/numFoam.txt');
order_T=1;
hold on
plt = plot(DisX(n+1,num_time+1-order_T*T_num:num_time-(order_T-1)*T_num),DisZ(n+1,num_time+1-order_T*T_num:num_time-(order_T-1)*T_num),'g-','LineWidth',2.5);
for j=num_time+1-order_T*T_num:DT:num_time-(order_T-1)*T_num%1:DT:num_time+1;%1:DT:T_num+1;%
plot(DisX(:,j),DisZ(:,j),'color',[39, 174, 96]/255,'Linewidth',0.1)
hold on
end
xlabel('x(m)')
ylabel('z(m)')
% legend(plt, 'Experiment data','Luhar & Nepf(2016)','Present Model','location','SouthEast','FontSize',8)
axis equal;
axis([-0.1 0.2 0 0.2])
grid on
figShape=[direction,'\Shape'];
saveas(gcf,figShape,'png');