function [A,B,C,D,Na_new,Nb_new,Nc_new,Nd_new]=Trimolecular_Reaction_Backward(A,B,C,D,Na,Nb,Nc,Nd,Pb,diff_2,diff_3,Da,Db,Dc,Lx)

warning('off','all')





r2=rand(Nd,1);
dtemp=(r2<Pb);
Dtemp=D(dtemp,:);
na=size(Dtemp,1);

temp=ones(na,3);
A=[A;temp];
B=[B;temp];
C=[C;temp];

for i=1:na
    
    xD=Dtemp(i,:);
    %Calculating position of new 3 molecules
    
    a1=pi*rand;
    a2=pi*rand;
    a3=pi*rand;
    a4=pi*rand;
    a5=2*pi*rand;
    
    x1=cos(a1)*diff_2;
    x2=sin(a1)*cos(a2)*diff_2;
    x3=sin(a1)*sin(a2)*cos(a3)*diff_2;
    x4=sin(a1)*sin(a2)*sin(a3)*cos(a4)*diff_3;
    x5=sin(a1)*sin(a2)*sin(a3)*sin(a4)*cos(a5)*diff_3;
    x6=sin(a1)*sin(a2)*sin(a3)*sin(a4)*sin(a5)*diff_3;
    
    eta1=xD';
    eta2=[x1;x2;x3];
    eta3=[x4;x5;x6];
    
    a=(1/Da+1/Db+1/Dc);
    b=(1/Da+1/Db);
    
    A1=[1/Da/a, 1/Db/a, 1/Dc/a;-1,1,0;-1/Da/b,-1/Db/b,1];
    xA=zeros(1,3);
    xB=zeros(1,3);
    xC=zeros(1,3);
    for j=1:3
        n=[eta1(j);eta2(j);eta3(j)];
        c=A1\n;
        xA(j)=c(1);
        xB(j)=c(2);
        xC(j)=c(3);
    end
    
    xA=mod(xA,Lx);
    xB=mod(xB,Lx);
    xC=mod(xC,Lx);
    
    
    
    A(Na+i,:)=xA;
    B(Nb+i,:)=xB;
    C(Nc+i,:)=xC;
end

D(dtemp,:)=[];


Na_new=size(A,1);
Nb_new=size(B,1);
Nc_new=size(C,1);
Nd_new=size(D,1);

