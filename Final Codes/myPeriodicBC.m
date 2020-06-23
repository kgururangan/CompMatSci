function [Kp,Mp] = myPeriodicBC(numnodes,Nx,Ny,K,M)
% K ane M must be in nine subvector set order!


% apply periodic b.c. (not adding bloch factors because u_k(r) is periodic
% in lattice (rve space)
% just mapping R -> L, T -> B, and 3 other corners -> Lb; interior nodes
% don't change
% number of dof associated with the left, bottom, and interior boundary
% nodes
nl = Ny - 1;    nb = Nx - 1;    ni = numnodes - 2*nl - 2*nb - 4;

nl = 2*nl;  nb = 2*nb;  ni = 2*ni;

            % bottom    % left     % LB       % int    
T_pbc = [   eye(nb,nb), zeros(nb,nl), zeros(nb,2), zeros(nb,ni);
            zeros(nl,nb),  eye(nl,nl), zeros(nl,2), zeros(nl,ni);
            zeros(2,nb), zeros(2,nl), eye(2,2), zeros(2,ni);
             zeros(ni,nb), zeros(ni,nl), zeros(ni,2), eye(ni,ni);
             % condensed out dof: ( T, R, RB, LT, RT)
             eye(nb,nb),  zeros(nb,nl), zeros(nb,2), zeros(nb,ni);
             zeros(nl,nb), eye(nl,nl), zeros(nl,2), zeros(nl,ni);
             zeros(2,nb), zeros(2,nl), eye(2,2), zeros(2,ni);
             zeros(2,nb), zeros(2,nl), eye(2,2), zeros(2,ni);
             zeros(2,nb), zeros(2,nl), eye(2,2), zeros(2,ni)];
         
 
Kp = T_pbc'*K*T_pbc;
Mp = T_pbc'*M*T_pbc;
             
             
    



end