function [A,B,Na_new,Nb_new]=Unimolecular_Reaction_Forward(A,B,Na,Pf)


r1=rand(1,Na);
dtemp=(r1<Pf);

if sum(dtemp)>0
    Btemp=A(dtemp,:);
    
%     n=size(Dtemp,1);
%     Btemp=zeros(n,3);
    
    %     parfor i=1:n
%     Btemp(1:n,:)=Dtemp(1:n,:);
    %     end
    B=[B;Btemp];
    
    A(dtemp,:)=[];
    
end
Na_new=size(A,1);
Nb_new=size(B,1);

