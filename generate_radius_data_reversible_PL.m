function yi = generate_radius_data_reversible_PL(slo,sinc,snum,N,alpha,PL)

warning('off','all')

steps = linspace(slo,slo + (snum-1)*sinc,snum);

error=1e-5*PL;
st = exp(steps);
if alpha<=1
    st=fliplr(st);
end
yi=zeros(1,length(st));
aa=ceil(alpha+1);
R = max(aa,2);


% options=optimset('Display','off');


Nr = 50*R+1;%number of nodes per unit length * R
r = linspace(0,R,Nr);
r=r';
sigma=alpha;

[~,r_alpha]=min(abs(r-sigma));
[~,r1] = min(abs(r-1));

gR1=max(r1,r_alpha);

dr = r(2)-r(1);
dt = 0.02*dr^2/2;
n = 3*(N-1);
t=0;
mu = (3*N - 5)/2;
surface = 2*(mu+1)*pi^(mu+1)/gamma(mu+2);

for scount = 1:length(st)
    count=1;
    if scount>1
        if yi(scount-1)>100 || isnan(yi(scount-1))
            yi(scount:length(st))=200;
            break
        end
    end
    s = st(scount);
    Dt = s^2/2;
    
    %calculate the K associated with the guess
    %construct Radial distribution function
    %     K=1./(r*r(r_alpha)).*1./(s*sqrt(2*pi)).*(exp(-(r-r(r_alpha)).^2/(2*s^2))-exp(-(r+r(r_alpha)).^2/(2*s^2)));
    %  K=r(r_alpha)./(surface*r*s*sqrt(2*pi)).*(exp(-(r-r(r_alpha)).^2/(2*s^2))-exp(-(r+r(r_alpha)).^2/(2*s^2)));
    % K=1/2*(1./r.^(mu)*r(r_alpha)^(1+mu)/s.*exp(-(r.^2+r(r_alpha)^2)/(4*s)).*besseli(-mu,r(r_alpha)*r/(2*s)));
    %     K=(1./r.^(mu)*r(r_alpha)^(1+mu)/s^2.*exp(-(r.^2+r(r_alpha)^2)/(2*s^2)).*besseli(mu,r(r_alpha)*r/(s)));
    %     K=(1./r.^(mu)*r(r_alpha)^(1+mu)/s.*exp(-(r.^2+r(r_alpha)^2)/(2*s^2)).*besseli(-mu,r(r_alpha)*r/(s)));
    %     K=1./(r*r(r_alpha)).*1./(s*sqrt(2*pi)).*(exp(-(r-r(r_alpha)).^2/(2*s^2))-exp(-(r+r(r_alpha)).^2/(2*s^2)));
    %          K_Reversible_Lipkova=r(r_alpha)./(r*s*sqrt(2*pi)).*(exp(-(r-r(r_alpha)).^2/(2*s^2))-exp(-(r+r(r_alpha)).^2/(2*s^2)));
    %          K_Reversible_Lipkova(1)=sqrt(2/pi)*r(r_alpha)^2/s^3*exp(-r(r_alpha)^2/(2*s^2));
    %          K=K_Reversible_Lipkova/r(r_alpha)^(2*mu+1);
    %          if alpha==0
    %              K_Reversible_Lipkova=sqrt(2/pi)/(4*pi*s^3)*exp(-r.^2/(2*s^2));
    %              K=K_Reversible_Lipkova;
    %          end
    
    %initial function
    amp=1/(2*s^2);
    delta=Dt*sqrt(amp/pi)*exp(-amp*(r-r(r_alpha)).^2)/r(r_alpha)^(2*mu+1);
%  K=1./((r*alpha).^mu).*exp(-(alpha-r).^2/(2*s^2)).*besseli(mu,r*alpha/s^2)*1/surface-1./((r*alpha).^mu).*exp(-(alpha+r).^2/(2*s^2)).*besseli(mu,r*alpha/s^2)*1/surface;
%  K(1)=1/(2)^mu*1/(Dt*gamma(1+mu))*exp(-alpha^2/(2*s^2))*(1/s^2)^mu*1/surface;
    if alpha==0
        r(r_alpha)=r(2);
% K=exp(-r.^2/(2*s^2))/(Dt*gamma(1+mu)*2^mu*s^(2*mu))*1/surface;
        delta=Dt*sqrt(amp/pi)*exp(-amp*(r-r(r_alpha)).^2)/r(r_alpha)^(2*mu+1);
    end
    g = ones(Nr,1);
        K=delta.*g;
   
%     K=0;
    
    reacted = 0;
    temp = PL*surface*(sum(r(1:r1).^(2*mu+1).*g(1:r1))-1/2*r(r1)^(2*mu+1)*g(r1))*dr; %The second half is the special stuff happening at r=1
    reactdiff = abs(reacted-temp)/temp;
    reacted = temp;
    w=reacted/(surface);
    
    g(gR1+1:end)=1;
% g(gR1+2:end)=1;
% g(gR1+1)=1/2*(g(gR1+2)+g(gR1));
    %     g=g+K*w;
    %     g(r_alpha+1:end)=1;
    
    %     g(1:r1-1) = zeros(r1-1,1);
    %     g(r_alpha+1:end)=1;
    g(1:r1-1)=(1-PL)*g(1:r1-1);
    
    
    g(r1)=(1-PL/2)*g(r1);
    
        g=g+w*K;
    
    
    temp = zeros(Nr,Nr);
    
    
    
    temp(Nr,Nr) = -2/dr;
    temp(Nr,Nr-1) = 1/dr - (n-1)/(2*R);
    for i = 2:Nr-1    %new i
        for j = 1:Nr    %old i
            
            if i == j
                
                temp(i,j) = -2/dr;
                
            elseif j == i-1
                
                
                temp(i,j) = 1/dr - (n-1)/(2*r(i));
                
            elseif j == i+1
                
                
                temp(i,j) = 1/dr + (n-1)/(2*r(i));
                
            end
            
        end
    end
    
    UpdateMatr = eye(Nr,Nr) + dt/dr*temp;
    
    
    
    while reactdiff > error
        while t < Dt
            
            g = UpdateMatr*g;
%             g(gR1+2:end)=1;
% g(gR1+1)=1/2*(g(gR1+2)+g(gR1));
          g(gR1+1:end)=1;
%           K=delta*dt;
%           g=g+w*K;
            g(1) = g(2);
            t = t+dt;
        end
          temp = PL*surface*(sum(r(1:r1).^(2*mu+1).*g(1:r1))-1/2*r(r1)^(2*mu+1)*g(r1))*dr;
%         temp = PL*surface*(sum(r(1:r1).^(2*mu+1).*g(1:r1))-1/2*r(r1)^(2*mu+1)*g(r1))*dr;
        reactdiff = abs(reacted-temp)/temp;
        reacted = temp;
        
        %         g=g+K*w;
        %         g(r_alpha+1:end)=1;
        %          g(r_alpha+2:end)=1;
        %                          g(r_alpha+1)=1/2*(g(r_alpha)+g(r_alpha+2));
        
        %                         g(r_alpha+1:end)=1;
        %          g(r1+1:end)=1;
%            K=delta.*g;
%   plot(r,g,'r')
%                         ylim([0 2])
%                         xlim([0 R])
%                         drawnow
%                         pause
        g(1:r1-1) = (1-PL)*g(1:r1-1);
        %      g(1:r1-1)=(1-PL)*g(1:r1-1);
        
        g(r1) = (1-PL/2)*g(r1);
        
        w=reacted/(surface);
      
       
        
                      
                         g=g+w*K;
        t=0;
        count=count+1;
        if isnan(reacted)
            yi(scount)=200;
            break
        end
        
        if alpha<1 && count>1e3
            yi(scount)=200;
            break
        end
        if count >1e4
            yi(scount)=200;
            break
        end
    end
%     count
%     reactdiff
    yi(scount) = reacted;
    if alpha<1 && count>1e3
        yi(scount)=200;
    end
      if isnan(reacted)
            yi(scount)=200;
      end
    
    
    fprintf('Alpha = %f,',alpha)
    fprintf(' for N = %f,', N)
    fprintf(' and steps = %f,', log(s))
    fprintf(' reduced reaction rate is %f \n', reacted)
end
if alpha<=1
    yi=fliplr(yi);
end
