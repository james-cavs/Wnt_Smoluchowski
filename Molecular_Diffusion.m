function A=Molecular_Diffusion(A,Na,move,L)
% warning('off','all')

%Start diffusion steps


if Na>0
    r1 = randn(Na,3); %random number to test SDE for each molecule
    %find the updated positions and the updated velocities
    A=A+move*r1;
    Lx = L(1);
    Ly = L(2);
    Lz = L(3);
    if (Lx == Ly) && (Lx == Lz)
        A=mod(A,Lx);
    else
        A(:,1)=mod(A(:,1),Lx);
        A(:,2)=mod(A(:,2),Ly);
        A(:,3)=mod(A(:,3),Lz);
    end
end