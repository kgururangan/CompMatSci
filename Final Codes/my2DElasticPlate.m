
clear all
close all


% Define side length of simulation box
L = 1;
% Periodic lattice spacing (distance between periodic structures)
a = [L,L];
% elastic properties of matrix and lattice material, respectively
youngs = [1,1];
poisson = [0,0];
% Need high density contrast (>100) to get "band structure" looking things
rho = [1,1];
% set mesh density
N = 100;
Nx = N;
Ny = N;
nen = 4;
numel = Nx*Ny;
numnodes = (Nx+1)*(Ny+1);
% generate mesh
elements = my2DMesh(L,Nx,Ny,youngs,poisson,rho);
hx = L/Nx;  hy = L/Ny;

Papp = 10;
%Pload = Papp*[1/Nx, 1/Ny];
Pload = Papp*[1/Nx, 1/Ny];

[K,M] = myKijMij_v2(elements,Nx,Ny);


% create the patch test force vector
F = zeros(2*numnodes,1);
idq = [];
tol = 1e-10;


for J = 1:numel
    
isboundary = elements{J,6};
X = elements{J,1};

% local to global coordinate mapping (1x4 vector)
sctr = [elements{J,2};
        elements{J,2}+numnodes];
% Implement mapping from element matrix ordering [u1x, u1y, u2x, u2y,..] to
% global matrix ordering [u1x, u2x, u3x, ... u1y, u2y, u3y,...]
%sctrVec = [sctr(1), sctr(1) + numnodes, sctr(2), sctr(2) + numnodes, sctr(3), sctr(3) + numnodes, sctr(4), sctr(4) + numnodes];
sctrVec = sctr(:)';

if isboundary
    
fe = zeros(2*nen,1);
    
% right hand side
if abs(X(2,1) - L) < tol
    fe(3) = Pload(1);  fe(5) = Pload(1);
    F(sctrVec) = fe;
    %idq = [idq, sctrVec(3), sctrVec(5)];
end

% top side
if abs(X(3,2) - L) < tol
    fe(6) = Pload(2);   fe(8) = Pload(2);
    F(sctrVec) = fe;
    %idq = [idq, sctrVec(6), sctrVec(8)];
end

% left side
%if X(1,1) == 0
%    fe(1) = -Pload(1);  fe(7) = -Pload(1);
%    F(sctrVec) = fe;
%    idq = [idq, sctrVec(1), sctrVec(7)];
%end
% bottom side
%if X(1,2) == 0
%    fe(2) = -Pload(2);  fe(4) = -Pload(2);
%    F(sctrVec) = fe;
%    idq = [idq, sctrVec(2), sctrVec(4)];
%end

               
end                
            
end

% constrain x and y of lower left corner & constrain y of lower right
% corner & constrain x of upper left corner
idd = [1, Ny*(Nx+1)+1, 1+numnodes, Nx+1+numnodes];

idf = 1:2*numnodes;
idf(find(idf==idd(1))) = [];
idf(find(idf==idd(2))) = [];
idf(find(idf==idd(3))) = [];
idf(find(idf==idd(4))) = [];

u = zeros(2*numnodes,1);

u(idf) = K(idf,idf)\( F(idf) - K(idf,idd)*u(idd) );
F(idd) = K(idd,idf)*u(idf) + K(idd,idd)*u(idd);

%% formatting and plotting displacements
close all

UX = reshape(u(1:numnodes),[Ny+1,Nx+1]);    UX = UX';
UY = reshape(u(numnodes+1:2*numnodes),[Ny+1,Nx+1]); UY = UY';

uxy = UX(:,floor(end/2));  uxx = UX(floor(end/2),:);
uyx = UY(floor(end/2),:);  uyy = UX(floor(end/2),:);

uxdiag = diag(UX);  uydiag = diag(UY);

ur = 1/sqrt(2)*(uxdiag + uydiag);
uth = 1/sqrt(2)*(-uxdiag + uydiag);

xe = 0:hx:L;
ye = 0:hy:L;

figure(1)
subplot(211)
imagesc(xe,ye,UX)
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
title('u_x','FontSize',14)
colorbar
subplot(212)
imagesc(xe,ye,UY)
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
title('u_y','FontSize',14)
colorbar

figure(2)
subplot(231)
plot(linspace(0,L,Nx+1),uxx,'Linewidth',2)
xlabel('x','FontSize',14)
title('u_x(x,y=L/2)','FontSize',14)
ll = legend('\epsilon_{xx} = 9.93');  set(ll,'FontSize',12,'Location','NorthWest');
grid on
subplot(232)
plot(linspace(0,L,Ny+1),uyy,'Linewidth',2)
xlabel('y','FontSize',14)
title('u_y(x=L/2,y)','FontSize',14)
ll = legend('\epsilon_{yy} = 9.93');  set(ll,'FontSize',12,'Location','NorthWest');
grid on
subplot(234)
plot(linspace(0,L,Ny+1),uxy,'Linewidth',2)
xlabel('y','FontSize',14)
title('u_x(x=L/2,y)','FontSize',14)
ll = legend('\epsilon_{xy} = 5.12');  set(ll,'FontSize',12,'Location','SouthEast');
grid on
subplot(235)
plot(linspace(0,L,Ny+1),uyx,'Linewidth',2)
xlabel('x','FontSize',14)
title('u_y(x,y=L/2)','FontSize',14)
ll = legend('\epsilon_{yx} = 5.12');  set(ll,'FontSize',12,'Location','SouthEast');
grid on
subplot(233)
plot(linspace(0,L*sqrt(2),Nx+1),ur,'Linewidth',2)
xlabel('r','FontSize',14)
title('u_r','FontSize',14)
ll = legend('\epsilon_{rr} = 13.79');  set(ll,'FontSize',12,'Location','SouthEast');
grid on
subplot(236)
plot(linspace(0,L*sqrt(2),Nx+1),uth,'Linewidth',2)
xlabel('r','FontSize',14)
title('u_{\theta}','FontSize',14)
ll = legend('\epsilon_{{\theta}r} = 0');  set(ll,'FontSize',12,'Location','SouthEast');
grid on


%% post-processing for strains - not working
% need B-matrix evaluated for every element (as a function of xi, eta ->
% x,y)
% interpolate strains over each element
% not so easy


EXX = zeros(Ny,Nx);     EYY = zeros(Ny,Nx);     EXY = zeros(Ny,Nx);
exx = zeros(numel,1);   eyy = exx;  exy = exx;
force_x = [];   force_y = [];
elcount = 1;

for i = Nx
    for j = 1:Ny

    
    conn = elements{elcount,2};
   
    %force_x = [force_x; F(conn)];
    %force_y = [force_y; F(conn+numnodes)];
    disp_x = u(conn);
    disp_y = u(conn+numnodes);
    
    d = [disp_x'; disp_y'];
    d = d(:);
    B = Bglob{elcount};
    
    strain = B*d;
    exx(elcount) = strain(1);
    eyy(elcount) = strain(2);
    exy(elcount) = 0.5*strain(3);  
    
    
    
    elcount = elcount + 1;
    
    EXX(j,i) = strain(1);
    EYY(j,i) = strain(2);
    EXY(j,i) = strain(3)*0.5;
    
    end
    
    
end
    

        
        
%EXX = reshape(exx,[Ny,Nx]);
%EYY = reshape(eyy,[Ny,Nx]);
%EXY = reshape(exy,[Ny,Nx]);

%FORCE_X = reshape(force_x,[Ny+1,Nx+1]);
%FORCE_Y = reshape(force_y,[Ny+1,Nx+1]);

%UX = reshape(u(1:numnodes),[Ny+1,Nx+1]);
%UY = reshape(u(numnodes+1:2*numnodes),[Ny+1,Nx+1]);

% element spatial dimension vector
Xe = 0:hx:L;
Ye = 0:hy:L;

exx_1 = EXX(floor(end/2),:);
eyy_1 = EYY(floor(end/2),:);
exy_1 = EXY(floor(end/2),:);

figure(1)
imagesc(Xe,Ye,EXX)
title('\epsilon_{xx}')
colormap('cool')
colorbar
figure(2)
imagesc(Xe,Ye,EYY)
title('\epsilon_{yy}')
colormap('cool')
colorbar
figure(3)
imagesc(Xe,Ye,EXY);
title('\epsilon_{xy}')
colormap('cool')
colorbar

figure(4)
subplot(3,1,1)
plot(exx_1);
xlabel('x')
title('\epsilon_{xx}')
grid on
subplot(3,1,2)
plot(eyy_1);
xlabel('x')
title('\epsilon_{yy}')
grid on
subplot(3,1,3)
plot(exy_1);
xlabel('x')
title('\epsilon_{xy}')
grid on
suptitle('Line through y = L/2')
    