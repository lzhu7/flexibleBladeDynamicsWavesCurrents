function J_value=J_value(Jacob,Eq_Govern_Done,Xi,Xj,u_wave_old,u_wave_new,v_wave_old,v_wave_new,n,m,w0,Umax)

%Boundary conditions
    J_value=zeros(6*(n+1),6*(n+1));
    J_value(1,3)=1;
    J_value(2,4)=1;
    J_value(3,6)=1;
    J_value(4,6*n+1)=1;
    J_value(5,6*n+2)=1;
    J_value(6,6*n+5)=1;
                
for k=2:n+1
    X0=Xi(6*(k-2)+1:6*k);
    X=Xj(6*(k-2)+1:6*k);
    Eq=@(X) Eq_Govern_Done(X0,X,u_wave_old(k-1:k),u_wave_new(k-1:k),v_wave_old(k-1:k),v_wave_new(k-1:k),m(k-1),w0(k-1),Umax(k-1));
    J_value(6*(k-1)+1:6*k,6*(k-2)+1:6*k)=Jacob(Eq,X);
end