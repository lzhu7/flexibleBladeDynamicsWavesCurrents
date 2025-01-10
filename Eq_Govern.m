function [Y,KC,Cdn,DragT,DragN,MaT,MaN,BT,BN,TT,TN,VBT,WT,AMN,VBN,WN]...
    =Eq_Govern(X0,X,Ux_wave_old,Ux_wave_new,Uz_wave_old,Uz_wave_new,m_unit,ma_unit,Dt,Ds,w0_unit,rhou_water,Cdt,Cdn,diameter,Area_cross,E,I,thick,T_wave,Umax)

%% govern equation
%g=9.81;
%X=[T,Q,u,v,Omega,phi,T,Q,u,v,Omega,phi]

%%
KC=T_wave/diameter*Umax;
if isnan(Cdn)
    Cdn=max(1.95,10*KC^(-1/3));     
    if KC<20
        Cm_1=(1+(0.35*KC^(2/3)));
    else
        Cm_1=(1+(0.15*KC^(2/3)));
    end
    Cm_2=(1+((KC-18)^2)/49);
    Cm=min(Cm_1,Cm_2);
    ma_unit=Cm*ma_unit;
end

Re=Umax*diameter/1e-6;
if isnan(Cdt)
    Cdt = 0.074 * Re^(-1/5);
end

%%
DragT=-1/2*rhou_water*Cdt*2*(diameter+thick)/4*( ...%drag tangent
        abs(X(9)-(Ux_wave_new(2)*cos(X(12))+Uz_wave_new(2)*sin(X(12))))*(X(9)-(Ux_wave_new(2)*cos(X(12))+Uz_wave_new(2)*sin(X(12)))) ...
        +abs(X(3)-(Ux_wave_new(1)*cos(X(6))+Uz_wave_new(1)*sin(X(6))))*(X(3)-(Ux_wave_new(1)*cos(X(6))+Uz_wave_new(1)*sin(X(6)))) ...
        +abs(X0(9)-(Ux_wave_old(2)*cos(X0(12))+Uz_wave_old(2)*sin(X0(12))))*(X0(9)-(Ux_wave_old(2)*cos(X0(12))+Uz_wave_old(2)*sin(X0(12)))) ...
        +abs(X0(3)-(Ux_wave_old(1)*cos(X0(6))+Uz_wave_old(1)*sin(X0(6))))*(X0(3)-(Ux_wave_old(1)*cos(X0(6))+Uz_wave_old(1)*sin(X0(6)))) ...
    );
MaT=m_unit/2/Dt*(X(9)+X(3)-X0(9)-X0(3))-m_unit/8/Dt*(X(10)+X(4)+X0(10)+X0(4))*(X(12)+X(6)-X0(12)-X0(6));
TT=-1/4*(X(11)*X(8)+X0(11)*X0(8)+X(5)*X(2)+X0(5)*X0(2));
BT=+1/2/Ds*(X(7)+X0(7)-X(1)-X0(1));
VBT=rhou_water*diameter*thick*(...
        1/2/Dt*(Ux_wave_new(2)+Ux_wave_new(1)-Ux_wave_old(2)-Ux_wave_old(1))...
        *1/4*(cos(X(12))+cos(X(6))+cos(X0(12))+cos(X0(6)))...
        +1/2/Dt*(Uz_wave_new(2)+Uz_wave_new(1)-Uz_wave_old(2)-Uz_wave_old(1))...
        *1/4*(sin(X(12))+sin(X(6))+sin(X0(12))+sin(X0(6)))...
    ); 
WT=-w0_unit/4*(sin(X(12))+sin(X(6))+sin(X0(12))+sin(X0(6)));

Y(1,1)=-MaT+BT+TT+DragT+VBT+WT;

%%
DragN=-1/2*rhou_water*Cdn*diameter/4*( ...%drag normal
        abs(X(10)-(-Ux_wave_new(2)*sin(X(12))+Uz_wave_new(2)*cos(X(12))))*(X(10)-(-Ux_wave_new(2)*sin(X(12))+Uz_wave_new(2)*cos(X(12)))) ...
        +abs(X(4)-(-Ux_wave_new(1)*sin(X(6))+Uz_wave_new(1)*cos(X(6))))*(X(4)-(-Ux_wave_new(1)*sin(X(6))+Uz_wave_new(1)*cos(X(6)))) ...
        +abs(X0(10)-(-Ux_wave_old(2)*sin(X0(12))+Uz_wave_old(2)*cos(X0(12))))*(X0(10)-(-Ux_wave_old(2)*sin(X0(12))+Uz_wave_old(2)*cos(X0(12)))) ...
        +abs(X0(4)-(-Ux_wave_old(1)*sin(X0(6))+Uz_wave_old(1)*cos(X0(6))))*(X0(4)-(-Ux_wave_old(1)*sin(X0(6))+Uz_wave_old(1)*cos(X0(6)))) ...
    );
MaN=m_unit/2/Dt*(X(10)+X(4)-X0(10)-X0(4))+m_unit/8/Dt*(X(9)+X(3)+X0(9)+X0(3))*(X(12)+X(6)-X0(12)-X0(6));
BN=+1/2/Ds*(X(8)+X0(8)-X(2)-X0(2));
TN=+1/4*(X(11)*X(7)+X0(11)*X0(7)+X(5)*X(1)+X0(5)*X0(1));
AMN=-ma_unit/2/Dt*(X(10)-(-Ux_wave_new(2)*sin(X(12))+Uz_wave_new(2)*cos(X(12)))...
    +X(4)-(-Ux_wave_new(1)*sin(X(6))+Uz_wave_new(1)*cos(X(6))) ...
    -(X0(10)-(-Ux_wave_old(2)*sin(X0(12))+Uz_wave_old(2)*cos(X0(12))))...
    -(X0(4)-(-Ux_wave_old(1)*sin(X0(6))+Uz_wave_old(1)*cos(X0(6))))); ...%added mass
    % added mass should be the derivative of the relative velocity, affect a lot
% AMNLuhar=-ma_unit/2/Dt*(X(10)+X(4)-X0(10)-X0(4))....
%     -ma_unit/8/Dt*(X(9)+X(3)+X0(9)+X0(3))*(X(12)+X(6)-X0(12)-X0(6))....  %%DV/Dt=pV/pt+wxV
%     -ma_unit/8/Dt*((Ux_wave_new(2)+Ux_wave_new(1)-Ux_wave_old(2)-Ux_wave_old(1))...
%     *(sin(X(12))+sin(X(6))+sin(X0(12))+sin(X0(6)))...
%     -(Uz_wave_new(2)+Uz_wave_new(1)-Uz_wave_old(2)-Uz_wave_old(1))....
%     *(cos(X(12))+cos(X(6))+cos(X0(12))+cos(X0(6))));...%%Luhar fomula
VBN=-rhou_water*diameter*thick*( ...
    -1/2/Dt*(Ux_wave_new(2)+Ux_wave_new(1)-Ux_wave_old(2)-Ux_wave_old(1))...
    *1/4*(sin(X(12))+sin(X(6))+sin(X0(12))+sin(X0(6)))...
    +1/2/Dt*(Uz_wave_new(2)+Uz_wave_new(1)-Uz_wave_old(2)-Uz_wave_old(1))...
    *1/4*(cos(X(12))+cos(X(6))+cos(X0(12))+cos(X0(6)))...
    );% ...%Virtual Buoyancy
WN=-w0_unit/4*(cos(X(12))+cos(X(6))+cos(X0(12))+cos(X0(6)));
Y(2,1)=-MaN+BN+TN+DragN+AMN+VBN+WN;

%%
EA = E*Area_cross;
Y(3,1)=-1/2/Dt*(X(7)+X(1)-X0(7)-X0(1))+EA/2/Ds*(X(9)+X0(9)-X(3)-X0(3))-EA/4*(X(11)*X(10)+X0(11)*X0(10)+X(5)*X(4)+X0(5)*X0(4));
%%
Y(4,1)=-1/2/Dt*(X(12)+X(6)-X0(12)-X0(6))+1/2/Ds*(X(10)+X0(10)-X(4)-X0(4))+1/4*(X(11)*X(9)+X0(11)*X0(9)+X(5)*X(3)+X0(5)*X0(3));
%%
Y(5,1)=E*I/Ds*(X(11)-X(5))+1/2*(X(8)+X(2));
%%
Y(6,1)=1/Ds*(X(12)-X(6))-1/2*(X(11)+X(5));