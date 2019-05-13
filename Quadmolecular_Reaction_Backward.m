function [A,B,C,D,E,Na_new,Nb_new,Nc_new,Nd_new,Ne_new]=Quadmolecular_Reaction_Backward(A,B,C,D,E,Na,Nb,Nc,Nd,Ne,Pb,diff_2,diff_3,diff_4,Da,Db,Dc,Dd,Lx)

warning('off','all')





r2=rand(1,Ne);
dtemp=(r2<Pb);
Etemp=E(:,dtemp);
na=size(Etemp,2);

temp=ones(3,na);
A=[A,temp];
B=[B,temp];
C=[C,temp];
D=[D,temp];

for i=1:na
    
    
    xE=Etemp(:,i);
    %Calculating position of new 3 molecules
    
    a1=pi*rand;
    a2=pi*rand;
    a3=pi*rand;
    a4=pi*rand;
    a5=pi*rand;
    a6=pi*rand;
    a7=pi*rand;
    a8=2*pi*rand;
    
    x1=cos(a1)*diff_2;
    x2=sin(a1)*cos(a2)*diff_2;
    x3=sin(a1)*sin(a2)*cos(a3)*diff_2;
    x4=sin(a1)*sin(a2)*sin(a3)*cos(a4)*diff_3;
    x5=sin(a1)*sin(a2)*sin(a3)*sin(a4)*cos(a5)*diff_3;
    x6=sin(a1)*sin(a2)*sin(a3)*sin(a4)*sin(a5)*cos(a6)*diff_3;
    x7=sin(a1)*sin(a2)*sin(a3)*sin(a4)*sin(a5)*sin(a6)*cos(a7)*diff_4;
    x8=sin(a1)*sin(a2)*sin(a3)*sin(a4)*sin(a5)*sin(a6)*sin(a7)*cos(a8)*diff_4;
    x9=sin(a1)*sin(a2)*sin(a3)*sin(a4)*sin(a5)*sin(a6)*sin(a7)*sin(a8)*diff_4;
    
    eta1=xE;
    eta2=[x1;x2;x3];
    eta3=[x4;x5;x6];
    eta4=[x7;x8;x9];
    
    a=(1/Da+1/Db+1/Dc+1/Dd);
    b=(1/Da+1/Db+1/Dc);
    c=(1/Da+1/Db);
    
    
    A1=[1/Da/a, 1/Db/a, 1/Dc/a, 1/Dd/a;-1,1,0,0;-1/Da/c,-1/Db/c,1,0;-1/Da/b,-1/Db/b,-1/Dc/b,1];
    xA=zeros(3,1);
    xB=zeros(3,1);
    xC=zeros(3,1);
    xD=zeros(3,1);
    for j=1:3
        n=[eta1(j);eta2(j);eta3(j);eta4(j)];
        c=A1\n;
        xA(j)=c(1);
        xB(j)=c(2);
        xC(j)=c(3);
        xD(j)=c(4);
    end
    
    
    
    xA=mod(xA,Lx);
    xB=mod(xB,Lx);
    xC=mod(xC,Lx);
    xD=mod(xD,Lx);
    
    A(:,Na+i)=xA;
    B(:,Nb+i)=xB;
    C(:,Nc+i)=xC;
    D(:,Nd+i)=xD;
    
    
    %     A=[A,xA];
    %     B=[B,xB];
    %     C=[C,xC];
    %     D=[D,xD];
    
    
    
end



    E(:,dtemp)=[];



Na_new=size(A,2);
Nb_new=size(B,2);
Nc_new=size(C,2);
Nd_new=size(D,2);
Ne_new=size(E,2);

