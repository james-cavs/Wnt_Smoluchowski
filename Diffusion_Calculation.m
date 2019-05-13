function [deltaN,DT]=Diffusion_Calculation(D,N)



Dtemp=0;
mtemp=0;
for i=1:N
    Dtemp=Dtemp+(1/D(i));
    if i==1
        continue
    end
    for m=1:i-1
        mtemp=mtemp+ 1/(D(i)*D(m));
    end
end
deltaN=Dtemp/mtemp;
DT=1;
Dt=zeros(1,N);
for i=2:N
    m=i-1;
    Dave=1/sum(1./D(1:m));
    Dt(i)=D(i)+Dave;
    DT=DT*(Dt(i)/deltaN)^(3/2);
end
