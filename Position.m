function    [DisX_new,DisZ_new,Length]=Position(DisX_old,DisZ_old,X0,X,Dt,Ds,n)

DisX_new=zeros(n+1,1);
DisZ_new=zeros(n+1,1);
for k=1:n+1
    DX=Dt*(1/2*(X(6*(k-1)+3)+X0(6*(k-1)+3)));%+1/Dt*(X(6*(k-1)+6)-X0(6*(k-1)+6))*Ds/2);%DisZ_old(k,1));
    DZ=Dt*(1/2*(X(6*(k-1)+4)+X0(6*(k-1)+4)));%+1/Dt*(X(6*(k-1)+6)-X0(6*(k-1)+6))*Ds/2);%DisX_old(k,1));
    DisX_new(k,1)=DisX_old(k,1)+DX.*cos(1/2*(X(6*(k-1)+6)+X0(6*(k-1)+6)))-DZ.*sin(1/2*(X(6*(k-1)+6)+X0(6*(k-1)+6)));
    DisZ_new(k,1)=DisZ_old(k,1)+DX.*sin(1/2*(X(6*(k-1)+6)+X0(6*(k-1)+6)))+DZ.*cos(1/2*(X(6*(k-1)+6)+X0(6*(k-1)+6)));
    
%     DisX_new(k,1)=DisX_old(k,1)+Dt*((X0(6*(k-1)+3).*cos(X0(6*(k-1)+6))-X0(6*(k-1)+4).*sin(X0(6*(k-1)+6)))+(X(6*(k-1)+3).*cos(X(6*(k-1)+6))-X(6*(k-1)+4).*sin(X(6*(k-1)+6))))/2;
%     DisZ_new(k,1)=DisZ_old(k,1)+Dt*((X0(6*(k-1)+3).*sin(X0(6*(k-1)+6))+X0(6*(k-1)+4).*cos(X0(6*(k-1)+6)))+(X(6*(k-1)+3).*sin(X(6*(k-1)+6))+X(6*(k-1)+4).*cos(X(6*(k-1)+6))))/2;
end

Length=0;
for k=2:n+1
    Length=Length+sqrt((DisX_new(k,1)-DisX_new(k-1,1))^2+(DisZ_new(k,1)-DisZ_new(k-1,1))^2);
end