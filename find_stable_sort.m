function [b,cbar_old]=find_stable_sort(x_count,x_index,box_count)


k=box_count;
c=zeros(k+1,1);
N=length(x_count);

for i=1:N
    c(x_count(i))=c(x_count(i))+1;
end

% temp = 1:(k+1);
% temp2 = x_count;

% c=histcounts(x_count,'BinLimits',[1 k]);
% if sum(c)>N
%     keyboard
% end



% c(end+1) = 0;






% keyboard
% keyboard
% parfor i = 1:(k+1)
%    c(i) = sum(temp2(:)==temp(i));
% end

% c(end)=0;
cbar=zeros(k+1,1);
cbar(1)=1;
for i=2:k+1
    cbar(i)=cbar(i-1)+c(i-1);
end
cbar(cbar>N)=N;
% cbar=cbar-1
cbar_old=cbar;
b=zeros(length(x_index),1);
% b_index=zeros(length(x_count),1);

% b(cbar(x_count(i))+1)=x_index(i);
for i=1:N
    %     b(cbar(x_count(i)))=x_count(i);
    b(cbar(x_count(i)))=x_index(i);
    %     b_index(i)=x_index(i);
    cbar(x_count(i))=cbar(x_count(i))+1;
end

% cbar(1:(end-1)) = cbar(2:end);
% b(1)=[];
% cbar(cbar==0)=1;
% b_index
% size(b_index)
% size(Atemp)
% keyboard

% A=Atemp(b_index,:);