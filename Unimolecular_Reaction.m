function [A,B,Na_new,Nb_new,Ntemp1,Ntemp2]=Unimolecular_Reaction(A,B,Na,Pf)


r1=rand(1,Na);
dtemp=(r1<Pf);
Dtemp=A(:,dtemp);
for i=1:size(Dtemp,2)
    B(:,end+1)=A(:,i);
end


A(:,dtemp)=[];

Na_new=size(A,2);
Nb_new=size(B,2);

Ntemp1=Na-Na_new;
Ntemp2=Nb_new-Nb;