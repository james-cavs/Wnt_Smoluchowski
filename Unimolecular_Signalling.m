function [A,Na_new]=Unimolecular_Signalling(Nb,A,Na_old,Pf)

r1=rand(1,Nb);
dtemp=(r1<Pf);

n=sum(dtemp);
% Atemp=zeros(n,3);
% for i=1:n
if n>0
    Atemp=rand(n,3);
    A=[A;Atemp];
end
% end

Na_new=size(A,1);
