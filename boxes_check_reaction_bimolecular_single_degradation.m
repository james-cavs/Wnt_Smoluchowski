function [A,B,Nb_new] = boxes_check_reaction_bimolecular_single_degradation(A,B,Pf,rho,box_count_x,box_count_y,box_count_z,spacing)

rho2 = rho^2;
Itemp=1;

ETA2=[];
Atemp=zeros(size(A,1),7);
Btemp=zeros(size(B,1),7);

Atemp(:,1:3)=A;
Btemp(:,1:3)=B;
% Atemp(4:6,:)=0;
% Btemp(4:6,:)=0;
Atemp(:,7)=1:size(A,1);
Btemp(:,7)=1:size(B,1);

% keyboard

Natemp=size(A,1);
Nbtemp=size(B,1);

% if Natemp <= Nbtemp
%     r1=rand(1,Natemp);
%     dtemp=(r1>=Pf);
%     Atemp(:,dtemp)=[];
%     Natemp=size(Atemp,2);
% elseif Natemp > Nbtemp
%     r1=rand(1,Nbtemp);
%     dtemp=(r1>=Pf);
%     Btemp(:,dtemp)=[];
%     Nbtemp=size(Btemp,2);
% end

temp=[Natemp,Nbtemp];
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
        
        
        
        
end



if ~isempty(Atemp)&&~isempty(Btemp)
    Atemp = boxes_update_allocation(Atemp,spacing);
    Btemp = boxes_update_allocation(Btemp,spacing);
    
    %  keyboard
    %
    %     [~,idx1]=sort(Atemp(:,4));
    %     [~,idx2]=sort(Btemp(:,4));
    %
    %     Atemp=Atemp(idx1,:);
    %     Btemp=Btemp(idx2,:);
    %
    %     [Atemp,~]=find_stable_sort(Atemp,box_count_x);
    
    [b_index,cbar]=find_stable_sort(Btemp(:,4),Btemp(:,7),box_count_x);
    
    Btemp=Btemp(b_index,:);
    
    % keyboard
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
                    Bitemp_start=cbar(a4(n)-1);
                    if cbar(a4(n))==cbar(a4(n)-1)
                        Bitemp_start=cbar(a4(n));
                        if cbar(a4(n))==cbar(a4(n)+1)
                            Bitemp_start=cbar(a4(n)+1);
                        end
                    end
                    
                end
                
                
                
                if x_count(2)~=box_count_x
                    Bitemp_end=cbar(a4(n)+2);
                    if cbar(a4(n)+2)==cbar(a4(n)+1)
                        Bitemp_end=cbar(a4(n)+1);
                        if cbar(a4(n))==cbar(a4(n)+1)
                            Bitemp_end=cbar(a4(n));
                        end
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
                for m=1:size(Btemp,1)
                    Px=(a1(n)-Btemp(m,1))^2;
                    Py=(a2(n)-Btemp(m,2))^2;
                    Pz=(a3(n)-Btemp(m,3))^2;
                    d1=Px+Py+Pz;
                    if (d1<=rho2)
                        ETA2=[ETA2;a7(n),Btemp(m,7)];
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
        
        if ~isempty(ETA3)
            
            if Itemp==2
                temp=ETA3(:,1);
                ETA3(:,1)=ETA3(:,2);
                ETA3(:,2)=temp;
            end
            
            
            
            
             B(ETA3(:,2),:)=[];
            
        end
    end
end



Nb_new=size(B,1);

