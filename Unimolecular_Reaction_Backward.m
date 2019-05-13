function [A,B,Na_new,Nb_new]=Unimolecular_Reaction_Backward(A,B,Na,Pf)

r1=rand(1,Na);
dtemp=(r1<Pf);
Dtemp=A(dtemp,:);

n=size(Dtemp,1);
Btemp=zeros(n,3);

for i=1:n
    Btemp(i,:)=A(i,:);
end
B=[B;Btemp];

A(dtemp,:)=[];

Na_new=size(A,1);
Nb_new=size(B,1);

