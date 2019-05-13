function [A,B,C,D,E,Na_new,Nb_new,Nc_new,Nd_new,Ne_new,Ntemp1,Ntemp2,Ntemp3,Ntemp4,Ntemp5] = boxes_check_reaction_quadmolecular(A,B,C,D,E,Pf,rho2temp,rho3temp,rho4temp,rho2,Da,Db,Dc,Dd,D2,D3,D4,deltaN,Ne,box_count_x,box_count_y,box_count_z,spacing)



Itemp=1;
ETA1x=[];
X2=[];
ETA2=[];
ETA3x=[];
ETA3x2=[];
ETA31=[];
ETA312=[];
ETA3=[];
ETA4=[];

Atemp=zeros(size(A,1),7);
Btemp=zeros(size(B,1),7);
Ctemp=zeros(size(C,1),7);
Dtemp=zeros(size(D,1),7);

Atemp(:,1:3)=A;
Btemp(:,1:3)=B;
Ctemp(:,1:3)=C;
Dtemp(:,1:3)=D;
% Atemp(4:6,:)=0;
% Btemp(4:6,:)=0;
Atemp(:,7)=1:size(A,1);
Btemp(:,7)=1:size(B,1);
Ctemp(:,7)=1:size(C,1);
Dtemp(:,7)=1:size(D,1);



Natemp=size(A,1);
Nbtemp=size(B,1);
Nctemp=size(C,1);
Ndtemp=size(D,1);

temp=[Natemp,Nbtemp,Nctemp,Ndtemp];
[~,I]=min(temp);

if I==2
    temp=Atemp;
    Atemp=Btemp;
    Btemp=temp;
    temp=Natemp;
    Natemp=Nbtemp;
    Nbtemp=temp;
    I=1;
    Itemp=2;
elseif I==3
    
    temp=Atemp;
    Atemp=Ctemp;
    Ctemp=temp;
    temp=Natemp;
    Natemp=Nctemp;
    Nctemp=temp;
    I=1;
    Itemp=3;
elseif I==4
    
    temp=Atemp;
    Atemp=Dtemp;
    Dtemp=temp;
    temp=Natemp;
    Natemp=Ndtemp;
    Ndtemp=temp;
    I=1;
    Itemp=4;
end


switch I
    case 1
        
        r1=rand(Natemp,1);
        dtemp=(r1>=Pf);
        Atemp(dtemp,:)=[];
        Natemp=size(Atemp,1);
    case 2
        
        r1=rand(Nbtemp,1);
        dtemp=(r1>=Pf);
        Btemp(dtemp,:)=[];
        Nbtemp=size(Btemp,1);
        
    case 3
        
        r1=rand(Nctemp,1);
        dtemp=(r1>=Pf);
        Ctemp(dtemp,:)=[];
        Nctemp=size(Ctemp,1);
        
    case 4
        
        r1=rand(Ndtemp,1);
        dtemp=(r1>=Pf);
        Dtemp(dtemp,:)=[];
        Ndtemp=size(Dtemp,1);
end




if ~isempty(Atemp)&&~isempty(Btemp)&&~isempty(Ctemp)&&~isempty(Dtemp)
    Atemp = boxes_update_allocation(Atemp,spacing);
    Btemp = boxes_update_allocation(Btemp,spacing);
    Ctemp = boxes_update_allocation(Ctemp,spacing);
    Dtemp = boxes_update_allocation(Dtemp,spacing);
    
    [~,idx1]=sort(Atemp(:,4));
    [~,idx2]=sort(Btemp(:,4));
    [~,idx3]=sort(Ctemp(:,4));
    [~,idx4]=sort(Dtemp(:,4));
    
    Atemp=Atemp(idx1,:);
    Btemp=Btemp(idx2,:);
    Ctemp=Ctemp(idx3,:);
    Dtemp=Dtemp(idx4,:);

    
    
    
    switch I
        case 1
            
            a1=Atemp(:,1);
            a2=Atemp(:,2);
            a3=Atemp(:,3);
            a4=Atemp(:,4);
            a5=Atemp(:,5);
            a6=Atemp(:,6);
            a7=Atemp(:,7);
            
            Btemptemp=Btemp;
            Ctemptemp=Ctemp;
            Dtemptemp=Dtemp;
            
            
            for n=1:Natemp
                Btemp=Btemptemp;
           
                
                x_count=(a4(n)-1):(a4(n)+1);
                y_count=(a5(n)-1):(a5(n)+1);
                z_count=(a6(n)-1):(a6(n)+1);
                
                if x_count(1)==0
                    x_count(1)=box_count_x;
                end
                
                if y_count(1)==0
                    y_count(1)=box_count_y;
                end
                
                if z_count(1)==0
                    z_count(1)=box_count_z;
                end
                
                if x_count(3)>box_count_x
                    x_count(3)=1;
                end
                
                if y_count(3)>box_count_y
                    y_count(3)=1;
                end
                
                if z_count(3)>box_count_z
                    z_count(3)=1;
                end
                
                if x_count(2)==box_count_x
                    Bitemp=(Btemp(:,4)==1);
                    temp=Btemp(Bitemp,:);
                    Btemp(Bitemp,:)=[];
                    Btemp=[Btemp;temp];
                    
                    Bitemp_end=find(Btemp(:,4)==1,1,'last');
                    if isempty(Bitemp_end)
                        Bitemp_end=find(Btemp(:,4)==x_count(2),1,'last');
                    end
                    
                    if isempty(Bitemp_end)
                        Bitemp_end=find(Btemp(:,4)==x_count(1),1,'last');
                    end
                elseif x_count(2)==1
                    Bitemp=(Btemp(:,4)==box_count_x);
                    temp=Btemp(Bitemp,:);
                    Btemp(Bitemp,:)=[];
                    Btemp=[temp;Btemp];
                    
                    Bitemp_start=find(Btemp(:,4)==box_count_x,1);
                    if isempty(Bitemp_start)
                        Bitemp_start=find(Btemp(:,4)==x_count(2),1);
                    end
                    if isempty(Bitemp_start)
                        Bitemp_start=find(Btemp(:,4)==x_count(3),1);
                    end
                end
                if x_count(2)~=1
                    Bitemp_start=find(Btemp(:,4)==x_count(1),1);
                    if isempty(Bitemp_start)
                        Bitemp_start=find(Btemp(:,4)==x_count(2),1);
                    end
                    if isempty(Bitemp_start)
                        Bitemp_start=find(Btemp(:,4)==x_count(3),1);
                    end
                end
                
                
                
                if x_count(2)~=box_count_x
                    Bitemp_end=find(Btemp(:,4)==x_count(3),1,'last');
                    
                    if isempty(Bitemp_end)
                        Bitemp_end=find(Btemp(:,4)==x_count(2),1,'last');
                    end
                    
                    if isempty(Bitemp_end)
                        Bitemp_end=find(Btemp(:,4)==x_count(1),1,'last');
                    end
                    
                end
                Btemp=Btemp(Bitemp_start:Bitemp_end,:);
                
                
                if ~isempty(Btemp)
                    Bjtemp=(Btemp(:,5)==y_count(1))|(Btemp(:,5)==y_count(2))|(Btemp(:,5)==y_count(3));
                    Btemp=Btemp(Bjtemp,:);
                end
                if ~isempty(Btemp)
                    Bktemp=(Btemp(:,6)==z_count(1))|(Btemp(:,6)==z_count(2))|(Btemp(:,6)==z_count(3));
                    Btemp=Btemp(Bktemp,:);
                end
                
                for o=1:size(Btemp,1)
                    Ctemp=Ctemptemp;
                    
                    if x_count(1)==0
                        x_count(1)=box_count_x;
                    end
                    
                    if y_count(1)==0
                        y_count(1)=box_count_y;
                    end
                    
                    if z_count(1)==0
                        z_count(1)=box_count_z;
                    end
                    
                    if x_count(3)>box_count_x
                        x_count(3)=1;
                    end
                    
                    if y_count(3)>box_count_y
                        y_count(3)=1;
                    end
                    
                    if z_count(3)>box_count_z
                        z_count(3)=1;
                    end
                    
                    if x_count(2)==box_count_x
                        Citemp=(Ctemp(:,4)==1);
                        temp=Ctemp(Citemp,:);
                        Ctemp(Citemp,:)=[];
                        Ctemp=[Ctemp;temp];
                        
                        Citemp_end=find(Ctemp(:,4)==1,1,'last');
                        if isempty(Citemp_end)
                            Citemp_end=find(Ctemp(:,4)==x_count(2),1,'last');
                        end
                        
                        if isempty(Citemp_end)
                            Citemp_end=find(Ctemp(:,4)==x_count(1),1,'last');
                        end
                    elseif x_count(2)==1
                        Citemp=(Ctemp(:,4)==box_count_x);
                        temp=Ctemp(Citemp,:);
                        Ctemp(Citemp,:)=[];
                        Ctemp=[temp;Ctemp];
                        
                        Citemp_start=find(Ctemp(:,4)==box_count_x,1);
                        if isempty(Citemp_start)
                            Citemp_start=find(Ctemp(:,4)==x_count(2),1);
                        end
                        if isempty(Citemp_start)
                            Citemp_start=find(Ctemp(:,4)==x_count(3),1);
                        end
                    end
                    if x_count(2)~=1
                        Citemp_start=find(Ctemp(:,4)==x_count(1),1);
                        if isempty(Citemp_start)
                            Citemp_start=find(Ctemp(:,4)==x_count(2),1);
                        end
                        if isempty(Citemp_start)
                            Citemp_start=find(Ctemp(:,4)==x_count(3),1);
                        end
                    end
                    
                    
                    
                    if x_count(2)~=box_count_x
                        Citemp_end=find(Ctemp(:,4)==x_count(3),1,'last');
                        
                        if isempty(Citemp_end)
                            Citemp_end=find(Ctemp(:,4)==x_count(2),1,'last');
                        end
                        
                        if isempty(Citemp_end)
                            Citemp_end=find(Ctemp(:,4)==x_count(1),1,'last');
                        end
                        
                    end
                    Ctemp=Ctemp(Citemp_start:Citemp_end,:);
                    
                    
                    if ~isempty(Ctemp)
                        Cjtemp=(Ctemp(:,5)==y_count(1))|(Ctemp(:,5)==y_count(2))|(Ctemp(:,5)==y_count(3));
                        Ctemp=Ctemp(Cjtemp,:);
                    end
                    if ~isempty(Ctemp)
                        Cktemp=(Ctemp(:,6)==z_count(1))|(Ctemp(:,6)==z_count(2))|(Ctemp(:,6)==z_count(3));
                        Ctemp=Ctemp(Cktemp,:);
                    end
                    
                    
                    for k=1:size(Ctemp,1)
                        Dtemp=Dtemptemp;
                        
                        if x_count(1)==0
                            x_count(1)=box_count_x;
                        end
                        
                        if y_count(1)==0
                            y_count(1)=box_count_y;
                        end
                        
                        if z_count(1)==0
                            z_count(1)=box_count_z;
                        end
                        
                        if x_count(3)>box_count_x
                            x_count(3)=1;
                        end
                        
                        if y_count(3)>box_count_y
                            y_count(3)=1;
                        end
                        
                        if z_count(3)>box_count_z
                            z_count(3)=1;
                        end
                        
                        if x_count(2)==box_count_x
                            Ditemp=(Dtemp(:,4)==1);
                            temp=Dtemp(Ditemp,:);
                            Dtemp(Ditemp,:)=[];
                            Dtemp=[Dtemp;temp];
                            
                            Ditemp_end=find(Dtemp(:,4)==1,1,'last');
                            if isempty(Ditemp_end)
                                Ditemp_end=find(Dtemp(:,4)==x_count(2),1,'last');
                            end
                            
                            if isempty(Ditemp_end)
                                Ditemp_end=find(Dtemp(:,4)==x_count(1),1,'last');
                            end
                        elseif x_count(2)==1
                            Ditemp=(Dtemp(:,4)==box_count_x);
                            temp=Dtemp(Ditemp,:);
                            Dtemp(Ditemp,:)=[];
                            Dtemp=[temp;Dtemp];
                            
                            Ditemp_start=find(Dtemp(:,4)==box_count_x,1);
                            if isempty(Ditemp_start)
                                Ditemp_start=find(Dtemp(:,4)==x_count(2),1);
                            end
                            if isempty(Ditemp_start)
                                Ditemp_start=find(Dtemp(:,4)==x_count(3),1);
                            end
                        end
                        if x_count(2)~=1
                            Ditemp_start=find(Dtemp(:,4)==x_count(1),1);
                            if isempty(Ditemp_start)
                                Ditemp_start=find(Dtemp(:,4)==x_count(2),1);
                            end
                            if isempty(Ditemp_start)
                                Ditemp_start=find(Dtemp(:,4)==x_count(3),1);
                            end
                        end
                        
                        
                        
                        if x_count(2)~=box_count_x
                            Ditemp_end=find(Dtemp(:,4)==x_count(3),1,'last');
                            
                            if isempty(Ditemp_end)
                                Ditemp_end=find(Dtemp(:,4)==x_count(2),1,'last');
                            end
                            
                            if isempty(Ditemp_end)
                                Ditemp_end=find(Dtemp(:,4)==x_count(1),1,'last');
                            end
                            
                        end
                        Dtemp=Dtemp(Ditemp_start:Ditemp_end,:);
                        
                        
                        if ~isempty(Dtemp)
                            Djtemp=(Dtemp(:,5)==y_count(1))|(Dtemp(:,5)==y_count(2))|(Dtemp(:,5)==y_count(3));
                            Dtemp=Dtemp(Djtemp,:);
                        end
                        if ~isempty(Dtemp)
                            Dktemp=(Dtemp(:,6)==z_count(1))|(Dtemp(:,6)==z_count(2))|(Dtemp(:,6)==z_count(3));
                            Dtemp=Dtemp(Dktemp,:);
                        end
                        
                        
                        
                        
                        %                 ntemp=min(size(Btemp,2),size(Dtemp,2));
                        
                        for m=1:size(Dtemp,1)
                            Px=(Atemp(n,1)-Btemp(o,1))^2;
                            Py=(Atemp(n,2)-Btemp(o,2))^2;
                            Pz=(Atemp(n,3)-Btemp(o,3))^2;
                            eta2=Px+Py+Pz;
                            Px=(Ctemp(k,1)-(Atemp(n,1)/Da+Btemp(o,1)/Db)/(1/Da+1/Db))^2;
                            Py=(Ctemp(k,2)-(Atemp(n,2)/Da+Btemp(o,2)/Db)/(1/Da+1/Db))^2;
                            Pz=(Ctemp(k,3)-(Atemp(n,3)/Da+Btemp(o,3)/Db)/(1/Da+1/Db))^2;
                            
                            eta3=Px+Py+Pz;
                            
                            Px=(Dtemp(m,1)-(Atemp(n,1)/Da+Btemp(o,1)/Db+Ctemp(k,1)/Dc)/(1/Da+1/Db+1/Dc))^2;
                            Py=(Dtemp(m,2)-(Atemp(n,2)/Da+Btemp(o,2)/Db+Ctemp(k,2)/Dc)/(1/Da+1/Db+1/Dc))^2;
                            Pz=(Dtemp(m,3)-(Atemp(n,3)/Da+Btemp(o,3)/Db+Ctemp(k,3)/Dc)/(1/Da+1/Db+1/Dc))^2;
                            
                            eta4=Px+Py+Pz;
                            
                            
                            d1=eta2*deltaN/D2+eta3*deltaN/D3+eta4*deltaN/D4;
                            
                            if (d1<=rho2)
                                ETA2=[ETA2;Atemp(n,7),Btemp(o,7),Ctemp(k,7),Dtemp(m,7)];
                            end
                        end
                    end
                    
                end
                
                
            end
    end
    
    
    
    
    
    
    if ~isempty(ETA2)
        ETA3=ETA2;
        
        
        ETA3=unique(ETA3,'rows');
        [~,ia,~]=unique(ETA3(:,1));
        ETA3=ETA3(ia,:);
        [~,ia,~]=unique(ETA3(:,2));
        ETA3=ETA3(ia,:);
        [~,ia,~]=unique(ETA3(:,3));
        ETA3=ETA3(ia,:);
        [~,ia,~]=unique(ETA3(:,4));
        ETA3=ETA3(ia,:);
        
        
        if ~isempty(ETA3)
            
              if Itemp==2
                temp=ETA3(:,1);
                ETA3(:,1)=ETA3(:,2);
                ETA3(:,2)=temp;
                
            elseif Itemp==3
                
                temp=ETA3(:,1);
                ETA3(:,1)=ETA3(:,3);
                ETA3(:,3)=temp;
            elseif Itemp==4
                
                temp=ETA3(:,1);
                ETA3(:,1)=ETA3(:,4);
                ETA3(:,4)=temp;
            end
            
            
            temp=ones(3,size(ETA3,1));
            E=[E;temp];
            for i=1:size(ETA3,1)
                
                E(Ne+i,:)=(A(ETA3(i,1),1:3)/Da+B(ETA3(i,2),1:3)/Db+C(ETA3(i,3),1:3)/Dc+D(ETA3(i,4),1:3)/Dc)/(1/(1/Da+1/Db+1/Dc+1/Dd));
                
                
                %Add a new molecule D and remove 1 of A,B,C
                
            end
            
            
            A(ETA3(:,1),:)=[];
            B(ETA3(:,2),:)=[];
            C(ETA3(:,3),:)=[];
            D(ETA3(:,4),:)=[];
            %             ETA3
        end
    end
end


Na_new=size(A,1);
Nb_new=size(B,1);
Nc_new=size(C,1);
Nd_new=size(D,1);
Ne_new=size(E,1);
Ntemp=size(ETA3,1);
Ntemp1=Ntemp;
Ntemp2=Ntemp;
Ntemp3=Ntemp;
Ntemp4=Ntemp;
Ntemp5=Ntemp;