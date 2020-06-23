clear all
clc
close all
% Define side length of simulation box
L = 0.02;
% Periodic lattice spacing (distance between periodic structures)
a = [L,L];
% elastic properties of matrix and lattice material, respectively
youngs = [10,1];
poisson = [0.34, 0.34];
% Need high density contrast (>100) to get "band structure" looking things
rho = [5,1];
% set mesh density
Nx = 9;
Ny = 9;
% generate mesh
elements = my2DMesh_small(L,Nx,Ny,youngs,poisson,rho);
% mesh{I,1} = 4x2 matrix of element I's nodal coordinates in set ordering
% down_left -> down_right -> up_right -> up_left (1->2->3->4)
% mesh{I,2} = element connectivity for local nodes in element I

% kspace mesh density
nkpt = 50;
    
% number of eigenfrequencies we will obtain
numnodes = (Nx+1)*(Ny+1);
ndof = 2*numnodes;

% resolution of k-space plotting grid
s = linspace(0,1,nkpt);

% high symmetry points in space
kG = [0,0]';
kX = [pi/a(1), 0]';
kM = [pi/a(1),pi/a(2)]';

% Define path vector in k-space
kGtoX = (kG + s.*(kX - kG));
k1 = norm(kGtoX(:,end)- kGtoX(:,1));
kXtoM = (kX + s.*(kM - kX));
k2 = norm(kXtoM(:,end) - kXtoM(:,1));
kMtoG = (kM + s.*(kG - kM));
k3 = norm(kMtoG(:,end)- kMtoG(:,1));

% plotting k-vectors, normalized to length in k-space
k1path = linspace(0,k1,nkpt);
k2path = linspace(k1,k1+k2,nkpt);
k3path = linspace(k1+k2,k1+k2+k3,nkpt);
ktot = k1+k2+k3;
dk = ktot/(3*nkpt);

kcoord = [k1path, k2path, k3path];
kvec = [kGtoX, kXtoM, kMtoG];

nl = Ny - 1;    nb = Nx - 1;    ni = numnodes - 2*nl - 2*nb - 4;
nl = 2*nl;  nb = 2*nb;  ni = 2*ni;
nfree = nl + nb + ni + 2;

%Kassemble = zeros(nfree,nfree,length(kcoord));
freq = zeros(nfree,length(kcoord));
orderVec = myNineSetOrder(numnodes,Nx,Ny);

tic
for i = 1:length(kcoord)
    
    [K,M] = myKMblochelement(kvec(:,i),elements,Nx,Ny);

    [Kp, Mp] = myPeriodicBC(numnodes,Nx,Ny,K(orderVec,orderVec),M(orderVec,orderVec));
    
    [V,D] = eigs(Kp,Mp,nfree);
    
    freq(:,i) = diag(D);
    
    %Kassemble(:,:,i) = Kp*Mp^-1;
   
    
end
toc


%[Vseq,freq] = eigenshuffle(Kassemble);


%%
close all

nmodes = nfree;

vm = sqrt(youngs(1)/rho(1));
omega = sqrt(abs(freq));
omega = omega*L./vm;

%omega(find(omega >= 4.1)) = omega(find(omega>=4.1)) + 1.2;

figure(1)
hold on
for j = 1:length(kcoord)
    
    %if j <= 60
    %    omega(find(omega(:,j) >= 6.51),j) = omega(find(omega(:,j) >= 6.51),j) + 0.8;
    %end
    %if j > 60 && j<= 80
    %    omega(find(omega(:,j) >= 7.31),j) = omega(find(omega(:,j) >= 7.31),j) + 0.8;
    %end
    %if j > 80 && j<= 120
    %    omega(find(omega(:,j) >= 7.151),j) = omega(find(omega(:,j) >= 7.151),j) + 0.8;
    %end
    plot(kcoord(j),omega(end-nmodes+1:end,j),'bs','MarkerSize',2)


end
hold off
xlabel('k','FontSize',20)
ylabel('v/v_{mat}','FontSize',20)
xticks([0 k1, k1+k2, k1+k2+k3])
xticklabels({'\Gamma','X','M','\Gamma'})
title('Transverse and Longitudinal','FontSize',20)
axis([0 k1+k2+k3 0 6])
grid on


