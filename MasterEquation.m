function [n,phi]=MasterEquation(k1,k2,n,N)



phi=zeros(1,length(n));
phi(1)=1;
n_max=max(n);
for i=1:n_max
    phi(i+1)=((k2/k1)^i)/(factorial(i)^(N-1))*phi(1)*nchoosek(n_max,i);
    i
end

sum1=sum(phi);
phi=phi/sum1;
% plot(n,phi,'.')