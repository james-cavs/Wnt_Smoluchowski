
clear all
close all
a=5;
sims=1;
t_end=2000;


dt=5*10^(-a);
% dt = 0.01;
t=0:dt:t_end;

r=1e-3;
V=4*pi*r^3/3;
V=V/1e6;
Lc=V^(1/3);
V=V*1e-6*6.022*1e23;
Lc=Lc*100;

betas0=0.92;
Ds0=1.46e-3;
BsDs0=1.86e-3;
D0=7.29e-4;
ApAx0=9.75e-4;
APC0=88.7;
Axin0=4.93e-4;
GSK0=50;
beta0=153;

Population=[APC0,Axin0,GSK0,beta0,Ds0,BsDs0,betas0,D0,ApAx0]*V;
Population=round(Population)

k=[0.182,1.82e-2,5e-2,.267,.133,9.09e-2,.909,50,120,206,206,.417,.423,2.57e-4,8.22e-5,.167,30,1200,1e-6,1e-4,10,0.01,50,0.01,0.01];
GSKc=50;
APCc=100;
TCFc=15;
DSHc=100;

k1=k(1);
k2=k(2);
k3=k(3);
k5=k(5);
k6=k(6);
k6b=k(7);
k7=k(22);
k7b=k(8)*k7;
k8=k(23);
k8b=k(9)*k8;
k9=k(10);
k12=k(13);
k14=k(15);



time_scale=1;

KbTri=k6b-k6b*k6*GSKc/(k7b+k6*GSKc);
KfTri=k6*k7/(k7b+k6*GSKc);
k8_modified=k8-k8*k8b/(k8b+k9);
k9_modified=k9*k8/(k8b+k9);
k8_modified=(k8_modified+k9_modified)/2;
off = 1;
if ~off
    W=1;
else
    W=0;
end

Kw=k3*k1*W*DSHc/(k1*W+k2);
K_Wnt=Kw-Kw*k6*GSKc/(k7b+k6*GSKc);

KfTri=KfTri/time_scale;
KbTri=KbTri/time_scale;
K_Wnt=K_Wnt/time_scale;
k8_modified=k8_modified/time_scale;
k8_modified=k8_modified/V;
KfTri=KfTri/(V^2)

k=k/time_scale


k(6)=k(6)/V;
k(22)=k(22)/V;
k(23)=k(23)/V;
k(24)=k(24)/V;
k(25)=k(25)/V;
k(15)=k(15)*V;
k(13)=k(13)*V;

k1=k(1);
k2=k(2);
k3=k(3);
k4=k(4);
k5=k(5);
k6=k(6);
k6b=k(7);
k7=k(22);
k7b=k(8)*k7;
k7=k(22);
k8=k(23);
k8b=k(9)*k8;
k8=k(23);
k9=k(10);
k10=k(11);
k11=k(12);
k12=k(13);
k13=k(14);
k14=k(15);
k15=k(16);
k16=k(24);
k16b=k(17)*k16;
k16=k(24);
k17=k(25);
k17b=k(18)*k17;
k17=k(25);



Di=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]*1e-8/Lc^2*60/time_scale;

N=length(Di);

Dap=Di(1);
Dax=Di(2);
Dg=Di(3);
Dapaxg=Di(4);
Dtc=Di(5);
Db=Di(6);
Dtb=Di(7);
Daxg=Di(8);
Dapax=Di(9);
Dbax=Di(10);
Dbap=Di(11);
Dapsaxsg=Di(12);
Dbapsaxsg=Di(13);
Dbsapsaxsg=Di(14);
Dbstar=Di(15);

D2 = Dax + Dap;
D3 = Dg + 1/(1/Dap+1/Dax);
N2=2;
N1=3;

Dspecial=Db+Dapsaxsg;
mu1=(3*N1-5/2);
mu2=(3*N2-5)/2;

Di1=[Dap,Dax,Dg];
Di2=[Dtc,Db];
[deltaN1,DT1]=Diffusion_Calculation(Di1,N1);
[deltaN2,DT2]=Diffusion_Calculation(Di2,N2);

P_lambda1=0.1;
P_lambda2=0.1;
P_lambda3=0.1;

lambda1=P_lambda1/dt;
lambda2=P_lambda2/dt;
lambda3=P_lambda3/dt;

alpha1=1.5;
alpha2=1.5;



slo1=-5;
sinc1=0.2;
s_end1=1;
snum1 = (s_end1-slo1)/sinc1+1;
% yi1 = generate_radius_data_reversible_PL(slo1,sinc1,snum1,N1,alpha1,P_lambda1);
yi1=[0.00344371048200078,0.00507029645345645,0.00717840063465215,0.0102378318494344,0.0151222463435550,0.0211575457061206,0.0303246505741743,0.0428438133011820,0.0597332700359851,0.0819422069153357,0.110748027237981,0.146767685598195,0.189379756835747,0.237450964134292,0.288274409313079,0.338288878655941,0.383979340533619,0.422847170679949,0.453790498360406,0.476995138864267,0.493441269184963,0.504418338525492,0.511218093317473,0.514977259804084,0.516671290793668,0.517195976933295,0.517282485035394,0.517287924419677,0.517288007780113,0.517288007939729,0.517288007939729];

slo2=-5;
sinc=0.2;
s_end2=1;
snum2 = (s_end2-slo2)/sinc+1;

% yi3=generate_radius_data_reversible_PL_Bi(slo2,sinc,snum2,N2,alpha2,P_lambda2);
yi3=[0.000887095202098135,0.00130582876754790,0.00185653700470257,0.00266891230822857,0.00399117333134329,0.00566341069393076,0.00827647342039037,0.0119784153308259,0.0172060497278603,0.0244844169446454,0.0346292460611835,0.0485283586074210,0.0669740709172454,0.0909611606068221,0.121058531113919,0.157119861889280,0.197872773587414,0.240818524351218,0.282607335574405,0.320007427408950,0.350873104357273,0.374531308409022,0.391547149908122,0.403167716858839,0.410726517847418,0.415317666150266,0.417745184533390,0.418708803549998,0.418937824632727,0.418962019064130,0.418962791952824];
% yi2=generate_radius_data_reversible_PL_Bi(slo2,sinc,snum2,N2,alpha2,P_lambda2);
yi2=yi3;


% APC +Axin + GSK x APC*/Axin*/GSK
rho1=find_radius(N1,alpha1,dt,DT1,deltaN1,KfTri,yi1,P_lambda1,slo1);

%APC*/Axin*/GSK + beta x APC*/Axin*/GSK/beta
rho3=find_radius(N2,alpha2,dt,DT2,deltaN2,k8_modified,yi3,P_lambda3,slo2);


rho=[rho1,rho3]
sigma=alpha2*rho;
sigma(1)=alpha1*rho1;
sigma
sigma1=sigma(1);
sigma3=sigma(2);


Pf1=lambda1*dt;
Pb1=1-exp(-(KbTri+K_Wnt)*dt);

Pf3=lambda3*dt;

P14=1-exp(-k14*dt);
P14=poissrnd(P14,1,length(t));
P15=1-exp(-k15*dt);
P4=1-exp(-k4*dt);
P5=1-exp(-k5*dt);
P12=k12*dt;
P12=poissrnd(P12,1,length(t));
P13=1-exp(-k13*dt);
P10=1-exp(-k10*dt);
P11=1-exp(-k11*dt);
P11=0;

NP=Population;



L=1;
X=[0,L;0,L;0,L];


x_min=X(1,1);
x_max=X(1,2);
y_min=X(2,1);
y_max=X(2,2);
z_min=X(3,1);
z_max=X(3,2);

Lx=x_max-x_min;
Ly=y_max-y_min;
Lz=z_max-z_min;

L=[Lx,Ly,Lz];
move=sqrt(2*Di(1)*dt);

rho12=rho(1)^2;
rho32=rho(2)^2;

rho12temp=(D2*rho12/deltaN1);
rho12t=sqrt(rho12temp);
rho13temp=rho12/deltaN1*D3;
rho13t=sqrt(rho13temp);

Nap_Out=zeros(1,sims);
Nax_Out=zeros(1,sims);
Ng_Out=zeros(1,sims);
Nb_Out=zeros(1,sims);
Nds_Out=zeros(1,sims);
Nbsds_Out=zeros(1,sims);
Nbs_Out=zeros(1,sims);
Nd_Out=zeros(1,sims);


n=length(t);

[~,~,~,box_length1,box_count_x1,box_count_y1,box_count_z1]=boxes_creation(Lx,Ly,Lz,rho1);
[~,~,~,box_length3,box_count_x3,box_count_y3,box_count_z3]=boxes_creation(Lx,Ly,Lz,rho3);

date=datestr(now,'yyyy_mm_dd_HH_MM_SS');
filename=sprintf('WntTriSimulation_Incomplete_%s.mat',date);
counter=round(length(t)/10);
for j=1:sims
     Nap=zeros(1,length(t));
    Nax=zeros(1,length(t));
    Ng=zeros(1,length(t));
    Nb=zeros(1,length(t));
    Nds=zeros(1,length(t));
    Nbsds=zeros(1,length(t));
    Nbs=zeros(1,length(t));
    Nd=zeros(1,length(t));
  
    
    Np=Population;
    Nap(1:2)=Np(1);
    Nax(1:2)=Np(2);
    Ng(1:2)=Np(3);
    Nb(1:2)=Np(4);
    Nds(1:2)=Np(5);
    Nbsds(1:2)=Np(6);
    Nbs(1:2)=Np(7);
    Nd(1:2)=Np(8);
     
    
    Ap = repmat(L,Nap(1),1).*rand(Nap(1),3).*ones(Nap(1),3);
    Ax = repmat(L,Nax(1),1).*rand(Nax(1),3).*ones(Nax(1),3);
    G = repmat(L,Ng(1),1).*rand(Ng(1),3).*ones(Ng(1),3);
    D = repmat(L,Nd(1),1).*rand(Nd(1),3).*ones(Nd(1),3);
    Ds = repmat(L,Nds(1),1).*rand(Nds(1),3).*ones(Nds(1),3);
    B = repmat(L,Nb(1),1).*rand(Nb(1),3).*ones(Nb(1),3);
    BsDs = repmat(L,Nbsds(1),1).*rand(Nbsds(1),3).*ones(Nbsds(1),3);
    diff_2=1/sqrt(deltaN1)*sqrt(D2)*sigma1;
    diff_3=1/sqrt(deltaN1)*sqrt(D3)*sigma1;
    

    td=0;


    for i=2:n
        Nap(i) = Nap(i-1);
         Nax(i) = Nax(i-1);
         Ng(i) = Ng(i-1);
                  Nd(i) = Nd(i-1);
         Nds(i) = Nds(i-1);
         Nb(i) = Nb(i-1);
         Nbsds(i) = Nbsds(i-1);
         Nbs(i) = Nbs(i-1);
         
        td=td+dt;
        
        %         Diffusion Step
        
        Ap = Molecular_Diffusion(Ap,Nap(i),move,L);
        Ax = Molecular_Diffusion(Ax,Nax(i),move,L);
        G = Molecular_Diffusion(G,Ng(i),move,L);
        
        D = Molecular_Diffusion(D,Nd(i),move,L);
        Ds = Molecular_Diffusion(Ds,Nds(i),move,L);
        B = Molecular_Diffusion(B,Nb(i),move,L);
        BsDs = Molecular_Diffusion(BsDs,Nbsds(i),move,L);

        
% %         Trimolecular reactions APC+Axin+GSK <> APC/Axin/GSK
        [Ap,Ax,G,D,Nap(i),Nax(i),Ng(i),Nd(i)]=boxes_check_reaction_trimolecular(Ap,Ax,G,D,Pf1,rho12temp,rho13temp,rho12,Dap,Dax,Dg,D2,D3,deltaN1,Nd(i-1),box_count_x1,box_count_y1,box_count_z1,box_length1);
        [Ap,Ax,G,D,Nap(i),Nax(i),Ng(i),Nd(i)]=Trimolecular_Reaction_Backward(Ap,Ax,G,D,Nap(i),Nax(i),Ng(i),Nd(i),Pb1,diff_2,diff_3,Dap,Dax,Dg,Lx);
       
        
        
% %         Bimolecular reaction APC*/Axin*/GSK+beta <>beta/APC*/Axin*/GSK


        [Ds,B,BsDs,Nds(i),Nb(i),Nbsds(i)]=boxes_check_reaction_bimolecular(Ds,B,BsDs,Pf3,rho3,Dapsaxsg,Db,Nbsds(i),box_count_x3,box_count_y3,box_count_z3,box_length3);
                
     %         Unimolecular reaction APC/Axin/GSK <> APC*/Axin*/GSK
        [D,Ds,Nd(i),Nds(i)]=Unimolecular_Reaction_Forward(D,Ds,Nd(i),P4);
        [Ds,D,Nds(i),Nd(i)]=Unimolecular_Reaction_Backward(Ds,D,Nds(i),P5);
        
        
        %         Reaction beta*/APC*/Axin*/GSK ->APC*/Axin*/GSK + beta*
        
        [BsDs,Ds,Nbsds(i),Nds(i),Nbs(i)]=Unimolecular_to_Bimolecular_Bs(BsDs,Ds,Nbs(i),Nbsds(i),P10);
        
        
        %         Unimolecular Beta decay
        [B,Nb(i)]=Unimolecular_Decay(B,Nb(i),P13);
        
                
        %         Unimolecular Axin decay
        [Ax,Nax(i)]=Unimolecular_Decay(Ax,Nax(i),P15);
        
        %         Zeroth order beta
        [B,Nb(i)]=Zeroth_Order(B,P12(i));
        
        %         Zeroth order Axin
        [Ax,Nax(i)]=Zeroth_Order(Ax,P14(i));

        Ndt=Nds(i)+Nbsds(i);
        Nbt=Nb(i)+Nbs(i)+Nbsds(i);
                Ndd1=Nbs(i-1);
                Ndd2=Nbs(i);
        if Ndd2>Ndd1
            td/time_scale
          
            SteadyState=[Nb(i),Nax(i),Ndt,Nbs(i)]
        end
        if mod(i,counter)==0
                    save(filename);
                    meow='meow'
        end
                
        
    end
    sims
    Nap_Out(j)=Nap(end);
    Nax_Out(j)=Nax(end);
    Ng_Out(j)=Ng(end);
        Nb_Out(j)=Nb(end);
        Nds_Out(j)=Nds(end);
       Nbsds_Out(j)=Nbsds(end);
    Nbs_Out(j)=Nbs(end);
    Nd_Out(j)=Nd(end);
    
    SteadyState=[Nb_Out(j),Nbs_Out(j),Nap_Out(j),Nax_Out(j),Ng_Out(j),Nd_Out(j)]
    
    
end

figure
plot(t/time_scale,Nd)
title('DC')
xlabel('Time (min)')
ylabel('Copy Number')

figure
plot(t/time_scale,Nb)
title('\beta-catenin')
xlabel('Time (min)')
ylabel('Copy Number')

date=datestr(now,'yyyy_mm_dd_HH_MM_SS');
filename=sprintf('WntTriSimulation_%s.mat',date);
save(filename);
