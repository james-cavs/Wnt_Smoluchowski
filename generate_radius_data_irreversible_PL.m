function yi = generate_radius_data_irreversible_PL(slo,sinc,snum,N,PL)

warning('off','all')

steps = linspace(slo,slo + (snum-1)*sinc,snum);


st = exp(steps);

yi=zeros(1,length(st));
R = 10;


% options=optimset('Display','off');


Nr = 20*R+1;%number of nodes per unit length * R
Nra = round(9*Nr/10);
r = linspace(0,R,Nr);
r=r';
% sigma=alpha;

% [~,r_alpha]=min(abs(r-sigma));
[~,r1] = min(abs(r-1));

gR1=max(r1);

dr = r(2)-r(1);
dt = 0.02*dr^2/2;
n = 3*(N-1);
t=0;
mu = (3*N - 5)/2;
surface = 2*(mu+1)*pi^(mu+1)/gamma(mu+2);
temp = zeros(Nr,Nr);

r_temp = r(Nra:Nr);
r0_temp = r(1:r1);
r_sum = sum(1./r_temp.^(4*mu));
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


parfor scount = 1:length(st)
    count=1;
    %     if scount>1
    %         if yi(scount-1)>100 || isnan(yi(scount-1))
    %             yi(scount:length(st))=200;
    %             break
    %         end
    %     end
    s = st(scount);
    Dt = s^2/2;
    
    
    g = ones(Nr,1);
    
    reacted = 0;
    temp = PL*surface*(sum(r0_temp.^(2*mu+1).*g(1:r1))-1/2*r(r1)^(2*mu+1)*g(r1))*dr; %The second half is the special stuff happening at r=1
    reactdiff = abs(reacted-temp)/temp;
    reacted = temp;
    
    g(1:r1-1)=(1-PL)*g(1:r1-1);
    g(r1)=(1-PL/2)*g(r1);
    
    
    
    
    a = -sum((1-g(Nra:Nr))./r_temp.^(2*mu))/sum(1./r_temp.^(4*mu));
    UpdateAdd = zeros(Nr,1);
    UpdateAdd(Nr) = dt/dr*(1/dr + (n-1)/(2*R))*(1+a/(R+dr)^(2*mu));
    
    
    while reactdiff > 0.00001
        t=0;
        

        while t < Dt
            g_temp = g(Nra:Nr);
            a = -sum((1-g_temp)./r_temp.^(2*mu))/r_sum;
            UpdateAdd(Nr) = dt/dr*(1/dr + (n-1)/(2*R))*(1+a/(R+dr)^(2*mu));
            g = UpdateMatr*g;
            g = g + UpdateAdd;
            g(1) = g(2);
            t = t+dt;
        end
        
        temp = PL*surface*(sum(r0_temp.^(2*mu+1).*g(1:r1))-1/2*r(r1)^(2*mu+1)*g(r1))*dr;
        reactdiff = abs(reacted-temp)/temp;
        reacted = temp;
        
        
        g(1:r1-1) = (1-PL)*g(1:r1-1);
        
        
        g(r1) = (1-PL/2)*g(r1);
        
        
        
        
        count=count+1;
        if isnan(reacted)
            yi(scount)=200;
            break
        end
        
        %         if alpha<1 && count>1e3
        %             yi(scount)=200;
        %             break
        %         end
        if count >1e4
            yi(scount)=200;
            break
        end
    end
    %     count
    %     reactdiff
    yi(scount) = reacted;
    %     if alpha<1 && count>1e3
    %         yi(scount)=200;
    %     end
    if isnan(reacted)
        yi(scount)=200;
    end
    
    
    %     fprintf('Alpha = %f,',alpha)
    fprintf(' for N = %f,', N)
    fprintf(' and steps = %f,', log(s))
    fprintf(' reduced reaction rate is %f \n', reacted)
end
% if alpha<=1
%     yi=fliplr(yi);
% end
