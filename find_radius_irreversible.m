function rho_smoldyn=find_radius_irreversible(N,dt,DT,deltaN,k1,yi,P_lambda,slo)
lambda=P_lambda/dt;

step=sqrt(2*dt*deltaN);
lo=0;
a=step;


while numrxnrate(step,a,N,DT,yi,lambda,deltaN,dt,slo) < k1*dt
   
    lo = a;
    a = 2*a;
    
end
dif = a - lo;

for n = 1:30
   
    dif = dif/2;
    a = lo + dif;
    if numrxnrate(step,a,N,DT,yi,lambda,deltaN,dt,slo) < k1*dt
        lo = a;
    end
    
end
rho_smoldyn = lo + 0.5*dif;
    
function k_guess = numrxnrate(step,a,N,DT,yi,lambda,deltaN,dt,slo)


P_lambda=lambda*dt;
mu = (3*N-5)/2;

% slo=-3;
sinc=0.2;
snum = length(yi);
step=log(step/a);
x=zeros(4,1);
z=zeros(4,1);

sindx = floor((step-slo)/sinc);
for i = 1:4
    x(i) = slo + (sindx - 2 + i)*sinc;
end
z(1) = (step-x(2))*(step-x(3))*(step-x(4))/(-6.0*sinc*sinc*sinc);
z(2) = (step-x(1))*(step-x(3))*(step-x(4))/(2.0*sinc*sinc*sinc);
z(3) = (step-x(1))*(step-x(2))*(step-x(4))/(-2.0*sinc*sinc*sinc);
z(4) = (step-x(1))*(step-x(2))*(step-x(3))/(6.0*sinc*sinc*sinc);
    
    
temp = 0;

for i = 1:4
   
    if (sindx-2+i >=snum)
      temp = temp + z(i)*pi^(mu+1)/gamma(mu+2)*P_lambda;
      
    elseif (sindx-2+i < 0)
            beta=a*sqrt(lambda/deltaN);
        
            if P_lambda<1
                A=1-2*mu*besseli(mu,beta)/(0.5*beta*(besseli(mu-1,beta)+besseli(mu+1,beta))+mu*besseli(mu,beta));
                %         Kpartial=2*mu*(mu+1)*alpha^(2*mu)*(1-A)/(alpha^(2*mu)-1+A);
                Kpartial=2*pi^(mu+1)*(A)/(gamma(mu));
            elseif P_lambda>=1
                %                Kpartial=2*alpha^(2*mu)*mu*(mu+1)/(alpha^(2*mu)-1);
                Kpartial=2*pi^(mu+1)/(gamma(mu));
            end
        
        temp = temp + z(i)*exp(2*x(i))*Kpartial;
        
    else
        temp = temp + z(i)*yi(sindx-1+i);
        
    end
    
end

k_guess = temp * a^(2*mu+2) * DT ;




