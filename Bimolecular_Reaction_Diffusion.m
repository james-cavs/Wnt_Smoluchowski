function A=Bimolecular_Reaction_Diffusion(A,Na,move,Lx)
warning('off','all')

%Start diffusion steps
    
    
    
    
    if Na>0
        r1 = randn(Na,3); %random number to test SDE for each molecule
        %find the updated positions and the updated velocities
        A=A+move*r1;
        A=mod(A,Lx);
    end