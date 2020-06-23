function [MESH] = my2DMesh_small(L,Nx,Ny,youngs,poisson,rho)

% currently only coded for a 9 x 9 element mesh 

% MESH = numelx6 cell array with entries
% col 1: element nodal coordinates in reference config
% col 2: element location matrix for local to global mapping
% col 3: element young's modulus
% col 4: element poisson ratio
% col 5: element density
% col 6: 1 or 0 indicating whether it's a boundary element

youngs_matrix = youngs(2);      youngs_per = youngs(1);
poisson_matrix = poisson(2);    poisson_per = poisson(1);
rho_matrix = rho(2);            rho_per = rho(1);


numel = Nx*Ny;
hx = L/Nx;  hy = L/Ny;

MESH = cell(numel,6);
% col 1: nodal coords 4x2
% col 2: connectivity matrix


% local node numbering is counter clockwise 1, 2, 3, 4 
loccoords = [0, 0;
             hx, 0;
             hx, hy;
             0, hy];                
locglob = [1, 2, 2+(Nx+1), 1+(Nx+1)];

%numnodes = (Nx+1)*(Ny+1);



c = 1;
for i = 1:Nx
    for j = 1:Ny
        
        %Le = zeros(4,numnodes);
        MESH{c,1} = [ loccoords(:,1) + hx*mod(i-1,Nx), loccoords(:,2) + hy*mod(j-1,Ny)];
        LM = locglob +  (j-1)*(Nx+1);
        %Le(1,LM(1)) = 1;
        %Le(2,LM(2)) = 1;
        %Le(3,LM(3)) = 1;
        %Le(4,LM(4)) = 1;
        MESH{c,2} = LM;
        
        % if element does not lie along the boundary, fill it with periodic
        % material
        if i > 3 && i < 7 && j > 3 && j < 7
            MESH{c,3} = youngs_per;
            MESH{c,4} = poisson_per;
            MESH{c,5} = rho_per;
            MESH{c,6} = 0;
        else
        % otherwise fill it with matrix (matrix = boundary)
            MESH{c,3} = youngs_matrix;
            MESH{c,4} = poisson_matrix;
            MESH{c,5} = rho_matrix;
            MESH{c,6} = 1;
        end
            
        
        
        
        
        %if mod(c,2) == 0
        %MESH{c,3} = youngs_matrix;
        %MESH{c,4} = poisson_matrix;
        %MESH{c,5} = rho_matrix;
        %else
        %MESH{c,3} = youngs_per;
        %MESH{c,4} = poisson_per;
        %MESH{c,5} = rho_per;
        %end
        c = c+1;
    end
    locglob = locglob + ones(1,4);

end
    
    
end