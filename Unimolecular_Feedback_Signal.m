function [A,B,Nb_new]=Unimolecular_Feedback_Signal(A,B,Na,Pf)



r1=rand(1,Na);
dtemp=(r1<Pf);

if sum(dtemp)>0
Btemp=A(dtemp,:);

% n=size(Dtemp,1);
% Btemp=zeros(n,3);

% for i=1:n
%     Btemp(1:n,:)=Dtemp(1:n,:);
% end
B=[B;Btemp];

end

Nb_new=size(B,1);

