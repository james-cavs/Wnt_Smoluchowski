function [A,Na_new]=Zeroth_Order(A,r1)



if r1>0
    
    Atemp=rand(r1,3);
    A=[A;Atemp];
end

Na_new=size(A,1);
