clear all
clc
close all
% Define side length of simulation box
L = 0.02;
% Periodic lattice spacing (distance between periodic structures)
a = [L,L];
% elastic properties of matrix and lattice material, respectively
youngs = [0.5,1];
poisson = [0.34, 0.34];
rho = [1,1];
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
freqL = zeros(nfree,length(kcoord));
freqT = freqL;
freqtot = freqL;
orderVec = myNineSetOrder(numnodes,Nx,Ny);

tic
for i = 1:length(kcoord)
    
    [Ktot,Mtot] = myKMblochelement(kvec(:,i),elements,Nx,Ny);
    [KL,KT,M] = myKMbloch_long_trans(kvec(:,i),elements,Nx,Ny);

    [KLp, Mp] = myPeriodicBC(numnodes,Nx,Ny,KL(orderVec,orderVec),M(orderVec,orderVec));
    [KTp, Mp] = myPeriodicBC(numnodes,Nx,Ny,KT(orderVec,orderVec),M(orderVec,orderVec));
    [Kptot, Mptot] = myPeriodicBC(numnodes,Nx,Ny,Ktot(orderVec,orderVec),Mtot(orderVec,orderVec));
    
    [VL,DL] = eigs(KLp,Mp,nfree);
    [VT,DT] = eigs(KTp,Mp,nfree);
    [Vtot,Dtot] = eigs(Kptot,Mptot,nfree);
    
    freqL(:,i) = diag(DL);
    freqT(:,i) = diag(DT);
    freqtot(:,i) = diag(Dtot);
    
    %Kassemble(:,:,i) = Kp*Mp^-1;
   
    
end
toc


%[Vseq,freq] = eigenshuffle(Kassemble);


%% Plotting
close all

nmodes = nfree;

vm = sqrt(youngs(2)/rho(2));

omegaL = sqrt(abs(freqL));
omegaL = omegaL*L./vm;

omegaT = sqrt(abs(freqT));
omegaT = omegaT*L./vm;

omegaTot = sqrt(abs(freqtot));
omegaTot = omegaTot*L./vm;


figure(1)
subplot(121)
hold on
for j = 1:length(kcoord)

    plot(kcoord(j),omegaT(end-nmodes+1:end,j),'ro','MarkerSize',4,'color',[0.8,0,0])


end
hold off
xlabel('$k$','Interpreter','latex')
ylabel('$v/v_{mat}$','Interpreter','latex')
xticks([0 k1, k1+k2, k1+k2+k3])
xticklabels({'\Gamma','X','M','\Gamma'})
title('Transverse Modes')
axis([0 k1+k2+k3 0 10])
set(gca,'FontSize',20,'Linewidth',2)
grid on

subplot(122)
hold on
for j = 1:length(kcoord)

    plot(kcoord(j),omegaL(end-nmodes+1:end,j),'ro','MarkerSize',4,'color',[0.8,0,0])


end
hold off
xlabel('$k$','Interpreter','latex')
ylabel('$v/v_{mat}$','Interpreter','latex')
xticks([0 k1, k1+k2, k1+k2+k3])
xticklabels({'\Gamma','X','M','\Gamma'})
title('Longitudinal Modes')
axis([0 k1+k2+k3 0 10])
set(gca,'FontSize',20,'Linewidth',2)
grid on

%%
figure(2)
hold on
for j = 1:length(kcoord)

    plot(kcoord(j),omegaTot(end-nmodes+1:end,j),'ro','MarkerSize',2,'color',[0.8,0,0])


end
hold off
xlabel('$k$','Interpreter','latex')
ylabel('$v/v_{mat}$','Interpreter','latex')
xticks([0 k1, k1+k2, k1+k2+k3])
xticklabels({'\Gamma','X','M','\Gamma'})
title('Linear Combination')
axis([0 k1+k2+k3 0 10])
set(gca,'FontSize',20,'Linewidth',2)
grid on

%suptitle('v^2_{per}/v^2_{mat} = 1.8')

