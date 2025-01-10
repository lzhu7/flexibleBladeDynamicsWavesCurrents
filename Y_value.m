function [Y_value,KC,Cdn,DragT,DragN,MaT,MaN,BT,BN,TT,TN,VBT,WT,AMN,VBN,WN]=Y_value(Eq_Govern_Done,Xi,Xj,u_wave_old,u_wave_new,v_wave_old,v_wave_new,n,m,w0,zv,Umax)

Cdn=zeros(n,1);
KC=zeros(n,1);

Y_value=zeros(6*(n+1),1);
Y_value(1,1)=Xj(3);
Y_value(2,1)=Xj(4);
Y_value(3,1)=Xj(6)-pi/2;
Y_value(4,1)=Xj(6*n+1);
Y_value(5,1)=Xj(6*n+2);
Y_value(6,1)=Xj(6*n+5);

for k=2:n+1
    X0=Xi(6*(k-2)+1:6*k);
    X=Xj(6*(k-2)+1:6*k);
    [Y_value(6*(k-1)+1:6*k,1),KC(k-1,1),Cdn(k-1,1),DragT(k-1,1),DragN(k-1,1),MaT(k-1,1),MaN(k-1,1),BT(k-1,1),BN(k-1,1),TT(k-1,1),TN(k-1,1),VBT(k-1,1),WT(k-1,1),AMN(k-1,1),VBN(k-1,1),WN(k-1,1)]=Eq_Govern_Done(X0,X,u_wave_old(k-1:k),u_wave_new(k-1:k),v_wave_old(k-1:k),v_wave_new(k-1:k),m(k-1),w0(k-1),Umax(k-1));
end