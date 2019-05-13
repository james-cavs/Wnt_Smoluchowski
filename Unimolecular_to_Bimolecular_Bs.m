function [A,B,Na_new,Nb_new,Nc_new]=Unimolecular_to_Bimolecular_Bs(A,B,Nc,Na,Pf)


r1=rand(Na,1);
dtemp=(r1<Pf);
n=sum(dtemp);
Dtemp=A(dtemp,:);
Btemp=zeros(n,3);


for i=1:n
    Btemp(i,:)=Dtemp(i,:);
    
end
B=[B;Btemp];
A(dtemp,:)=[];

Na_new=size(A,1);
Nb_new=size(B,1);
Nc_new=Nc+n;
