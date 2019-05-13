function [n,phi]=MasterEquation_Multiple_PR(k1,k2,n,N,M)

% N = order of LHS reaction
% M = order of RHS reaction

phi=zeros(1,length(n));
phi(1)=1;
n_max=max(n);
for i=1:n_max
%     phi(i+1)=((k2/k1)^i)/(factorial(i)^(N-1))*phi(1)*nchoosek(n_max,i);
    phi(i+1)=((k2/k1)^i)*factorial(n_max)^M/(factorial(i)^N*factorial(n_max-i)^M)*phi(1);
%     i
end

sum1=sum(phi);
phi=phi/sum1;
% plot(n,phi,'.')