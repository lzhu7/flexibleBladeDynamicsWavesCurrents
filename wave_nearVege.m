function [u_wave,v_wave,Umax]=wave_nearVege(DisX,DisZ,time,k_wave_number,Omega_wave)

load('./LuharData/PIV_Vel_f05a40.mat','UU','VV','U','y','t','Tw');

% Tw=2;
% The different variables are:
% UU, VV:                The horizontal and vertical velocities, respectively, over a wave period [dimensionless]
% U:                     The velocity scale used to normalize UU, VV [cm/s]
% y:                     Vertical positions for the PIV  data [cm]
% t:                     Normalized time for the PIV data, i.e. it goes from 0 to 2 \pi [rad]
% f:                     Wave Frequency [Hz]
% Tw:                    Wave Period [s]
[n_s,~]=size(DisX);
n_s=n_s-1;

[n_y,n_t]=size(UU);
uDisZ=zeros(n_s+1,n_t);
vDisZ=zeros(n_s+1,n_t);
u_wave=zeros(n_s+1,1);
v_wave=zeros(n_s+1,1);
Umax=zeros(n_s+1,1);

z_X0=zeros(n_y+1,1);
uX0=zeros(n_y+1,n_t);
vX0=zeros(n_y+1,n_t);
uX0(2:end,:)=UU*U/100;
uX0(1,:)=uX0(2,:);
vX0(2:end,:)=VV*U/100;
vX0(1,:)=vX0(2,:);
z_X0(2:end)=y/100;
z_X0(1)=z_X0(2)+0.031; %as the measurements is less than 0.2
tWave=zeros(n_t,1);
tWave(2:n_t,1)=t(2:n_t)/2/pi*Tw;

tDisX=time-k_wave_number/Omega_wave*DisX;
% tDisX=time-k_wave_number/Omega_wave*zeros(n_s+1,1);%Luhar
tDisX=tDisX-floor(tDisX/Tw)*Tw;

for i=1:n_t
    uDisZ(:,i) = interp1(z_X0,uX0(:,i),DisZ);
    vDisZ(:,i) = interp1(z_X0,vX0(:,i),DisZ);
end

for j=1:n_s+1
    u_wave(j,1)=interp1(tWave,uDisZ(j,:),tDisX(j));
    v_wave(j,1)=interp1(tWave,vDisZ(j,:),tDisX(j));
    Umax(j,1)=max(abs(uDisZ(j,:)));
end

% Umax=max(max(abs(uDisZ)));
% Umax=max(abs(u_wave));
