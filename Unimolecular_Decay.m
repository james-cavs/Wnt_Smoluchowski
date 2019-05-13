function [A,Na_new]=Unimolecular_Decay(A,Na,Pf)


r1=rand(Na,1);
dtemp=(r1<Pf);

A(dtemp,:)=[];

Na_new=size(A,1);
