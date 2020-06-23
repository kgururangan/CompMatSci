function[K,M] = myKMblochelement(kvec,elements,Nx,Ny)

k1 = kvec(1);   k2 = kvec(2);

numnodes = (Nx+1)*(Ny+1);
numel = Nx*Ny;
nen = 4;

% shape functions and derivatives on isoparametric domain
n1 = @(xi) 1/4*(1-xi(1)).*(1-xi(2));    dn11 = @(xi) -1/4*(1-xi(2));   dn12 = @(xi) -1/4*(1-xi(1));  
n2 = @(xi) 1/4*(1+xi(1)).*(1-xi(2));    dn21 = @(xi) 1/4*(1-xi(2));    dn22 = @(xi) -1/4*(1+xi(1));
n3 = @(xi) 1/4*(1+xi(1)).*(1+xi(2));    dn31 = @(xi) 1/4*(1+xi(2));    dn32 = @(xi) 1/4*(1+xi(1));
n4 = @(xi) 1/4*(1-xi(1)).*(1+xi(2));    dn41 = @(xi) -1/4*(1+xi(2));   dn42 = @(xi) 1/4*(1-xi(1));

% gaussian quadrature points
xg = [-1/sqrt(3), -1/sqrt(3);
      1/sqrt(3), -1/sqrt(3);
      1/sqrt(3), 1/sqrt(3);
      -1/sqrt(3), 1/sqrt(3)];
wt = [1, 1, 1, 1];
ng = length(wt);


K = zeros(2*numnodes,2*numnodes);
M = zeros(2*numnodes,2*numnodes);

%K = sparse(2*numnodes,2*numnodes);
%M = sparse(2*numnodes,2*numnodes);

% loop over the elements and assemble global matrices
for jj = 1:numel
    
    % preallocate element matrices
    me = zeros(2*nen,2*nen);
    ke = zeros(2*nen,2*nen);
    
    % element coordinates in (X,Y) reference configuration (4x2 matrix)
    X = elements{jj,1};
    
    % element mechanical properties
    E = elements{jj,3};             % young's modulus
    v = elements{jj,4};             % poisson ratio
    rho = elements{jj,5};           % density
    % lame parameters
    lambda = (v*E)/( (1+v)*(1-2*v) );
    mu = E/(2*(1+v));
    
    % construct jacobian matrix to get from isoparametric to reference
    % confi
    jacobian = @(xig) [X(1,1)*dn11(xig) + X(2,1)*dn21(xig) + X(3,1)*dn31(xig) + X(4,1)*dn41(xig),   X(1,1)*dn12(xig) + X(2,1)*dn22(xig) + X(3,1)*dn32(xig) + X(4,1)*dn42(xig);
                       X(1,2)*dn11(xig) + X(2,2)*dn21(xig) + X(3,2)*dn31(xig) + X(4,2)*dn41(xig),   X(1,2)*dn12(xig) + X(2,2)*dn22(xig) + X(3,2)*dn32(xig) + X(4,2)*dn42(xig)];
                   
    % local to global coordinate mapping (1x4 vector)
    sctr = [elements{jj,2};
            elements{jj,2}+numnodes];
    % Implement mapping from element matrix ordering [u1x, u1y, u2x, u2y,..] to
    % global matrix ordering [u1x, u2x, u3x, ... u1y, u2y, u3y,...]
    % sctrVec = [sctr(1), sctr(1) + numnodes, sctr(2), sctr(2) + numnodes, sctr(3), sctr(3) + numnodes, sctr(4), sctr(4) + numnodes];
    % sctrVec gets elements in K, M, F that line up with Voigt notation
    % ordering in element vectors
    sctrVec = sctr(:)';
    
            % Loop over integration points
            for g = 1:ng
               
               % gauss integration point in isoparametric domain
               xig = xg(g,:);

               jac = jacobian(xig);
               jac = jac';
               jacinv = jac^-1;
               JX = det(jac);
               
               % use jacobian to express shape function derivatives over
               % reference domain in terms of isoparametric derivatives
               dn1x = dn11(xig)*jacinv(1,1) + dn12(xig)*jacinv(2,1);
               dn1y = dn11(xig)*jacinv(1,2) + dn12(xig)*jacinv(2,2);
               
               dn2x = dn21(xig)*jacinv(1,1) + dn22(xig)*jacinv(2,1);
               dn2y = dn21(xig)*jacinv(1,2) + dn22(xig)*jacinv(2,2);
               
               dn3x = dn31(xig)*jacinv(1,1) + dn32(xig)*jacinv(2,1);
               dn3y = dn31(xig)*jacinv(1,2) + dn32(xig)*jacinv(2,2);
               
               dn4x = dn41(xig)*jacinv(1,1) + dn42(xig)*jacinv(2,1);
               dn4y = dn41(xig)*jacinv(1,2) + dn42(xig)*jacinv(2,2);
               
               %dn1x = 2/(L/N)*(dn11(xig) + dn12(xig));
               %dn1y = 2/(L/N)*(dn11(xig) + dn12(xig));
               %dn2x = 2/(L/N)*(dn21(xig) + dn22(xig));
               %dn2y = 2/(L/N)*(dn21(xig) + dn22(xig));
               %dn3x = 2/(L/N)*(dn31(xig) + dn32(xig));
               %dn3y = 2/(L/N)*(dn31(xig) + dn32(xig));
               %dn4x = 2/(L/N)*(dn41(xig) + dn42(xig));
               %dn4y = 2/(L/N)*(dn41(xig) + dn42(xig));
               
               
               % Construct the B-matrix in Voigt notation [dn1/dx, 0; 0, dn1/dy, dn1/dy, dn1/dx] 
               % evaluate on isoparametric domain, use jacobian to change
               % integration coordinate
               Be1 = [dn1x, 0;
                      0, dn1y;
                      dn1y, dn1x];
                  
               Bek1 = [1i*k1*n1(xig), 0;
                      0,  1i*k2*n1(xig);
                      1i*k2*n1(xig), 1i*k1*n1(xig)];
                  
               Be2 = [dn2x, 0;
                      0, dn2y;
                      dn2y, dn2x];
                  
               Bek2 = [1i*k1*n2(xig), 0;
                        0, 1i*k2*n2(xig);
                        1i*k2*n2(xig), 1i*k1*n2(xig)];
                    
               Be3 = [dn3x, 0;
                      0, dn3y;
                      dn3y, dn3x];
                  
               Bek3 = [1i*k1*n3(xig), 0;
                      0, 1i*k2*n3(xig);
                      1i*k2*n3(xig), 1i*k1*n3(xig)];

               Be4 = [dn4x, 0;
                      0, dn4y;
                      dn4y, dn4x];
                  
               Bek4 = [1i*k1*n4(xig), 0;
                      0, 1i*k2*n4(xig);
                      1i*k2*n4(xig), 1i*k1*n4(xig)];
                    
               % element B-matrix    
               Be = [Be1 + Bek1, Be2 + Bek2, Be3 + Bek3, Be4 + Bek4];
               
               % element nodal matrix: orders in {u1x, u1y, u2x, u2y, ... }
               Ne = [ n1(xig), 0, n2(xig), 0, n3(xig), 0, n4(xig), 0;
                     0, n1(xig), 0, n2(xig), 0, n3(xig), 0, n4(xig)];
               
               % compute element mass matrix
               me = me + wt(g)*(Ne')*rho*Ne*JX;
               
               % compute element stiffness matrix
               %C = [lambda + 2*mu, lambda, 0;
               %         lambda, lambda+2*mu, 0;
               %         0,          0,       mu];
               % plane strain
               %C = lambda/v * [1-v, v, 0;
               %                 v, 1-v, 0;
               %                 0,  0, (1-2*v)/2];
               % plane stress
               C = E/(1-v^2)*[ 1 v 0; 
                               v 1 0;
                               0 0 (1-v)/2];
                           
               ke = ke + wt(g)*Be'*C*Be*JX;
      
                
            end
            
            % scatter into global matrices
            K(sctrVec,sctrVec) = K(sctrVec,sctrVec) + ke;
            M(sctrVec,sctrVec) = M(sctrVec,sctrVec) + me;
            
                
   
end

end