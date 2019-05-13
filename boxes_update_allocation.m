function A = boxes_update_allocation(A,h)

A(:,4)=ceil(A(:,1)/h);
A(:,5)=ceil(A(:,2)/h);
A(:,6)=ceil(A(:,3)/h);

