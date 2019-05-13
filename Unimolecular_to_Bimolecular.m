function [A,B,C,Na_new,Nb_new,Nc_new]=Unimolecular_to_Bimolecular(A,B,C,Na,Pf)


r1=rand(1,Na);
dtemp=(r1<Pf);
Dtemp=A(:,dtemp);

n=size(Dtemp,1);
Btemp=zeros(n,3);


for i=1:n
    Btemp(i,:)=Dtemp(i,:);
    
end
B=[B;Btemp];
A(dtemp,:)=[];

Na_new=size(A,2);
Nb_new=size(B,2);
Nc_new=Nc+n;
