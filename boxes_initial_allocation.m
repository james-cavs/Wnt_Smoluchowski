function A = boxes_initial_allocation(A,box_count_x,box_count_y,box_count_z,x_box,y_box,z_box)

N=size(A,2);

for i=1:N
    for j=1:box_count_x
        if A(1,i)>=x_box(j)&&A(1,i)<x_box(j+1)
            A(4,i)=j;
            break
        end
    end
    
    for j=1:box_count_y
        if A(2,i)>=y_box(j)&&A(2,i)<y_box(j+1)
            A(5,i)=j;
            break
        end
    end
    
    for j=1:box_count_z
        if A(3,i)>=z_box(j)&&A(3,i)<z_box(j+1)
            A(6,i)=j;
            break
        end
    end
end
