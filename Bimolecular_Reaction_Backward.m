function [A,B,C,Na_new,Nb_new,Nc_new]=Bimolecular_Reaction_Backward(A,B,C,Na,Nb,Nc,Pb,sigma,Da,Db,Lx)
warning('off','all')





Ntemp=0;

r2=rand(Nc,1);
dtemp=(r2<Pb);
Dtemp=C(dtemp,:);
na=size(Dtemp,1);
for i=1:na
    
    Ntemp=Ntemp+1;
    xD=Dtemp(i,:);
    %Calculating position of new 3 molecules
    
    a1=pi*rand;
    a2=2*pi*rand;
    
    
    x1=cos(a1)*sigma;
    x2=sin(a1)*cos(a2)*sigma;
    x3=sin(a1)*sin(a2)*sigma;
    
    
    eta1=xD';
    eta2=[x1;x2;x3];
    
    a=(1/Da+1/Db);
    
    A1=[1/Da/a, 1/Db/a;-1,1];
    xA=zeros(1,3);
    xB=zeros(1,3);
    for j=1:3
        n=[eta1(j);eta2(j)];
        c=A1\n;
        xA(j)=c(1);
        xB(j)=c(2);
    end
    
    xA=mod(xA,Lx);
    xB=mod(xB,Lx);
    
    
    
    A(Na+i,:)=xA;
    B(Nb+i,:)=xB;
    
end

C(dtemp,:)=[];



Na_new=size(A,1);
Nb_new=size(B,1);
Nc_new=size(C,1);
