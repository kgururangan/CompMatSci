%% Simulated Annealing Reconstruction (Temperature Scheme)
clear all
clc
close all
%%%%%%%%%%%%%%%%%%%% Simulated Annealing Parameters %%%%%%%%%%%%%%%%%%%%%%

% Define the initial temperature we want to start at and minimum
% temperature
% Starting temperature must be defined by T0 = dE/log(2) where dE is an
% initial energy difference between image and target
T = 7e-5;
Tmin = 1e-10;
Etol = 1e-6;

% Annealing schedule parameter
lambda = 0.999;

% Number of trials we want to take for each temperature
L = 100;

% Initialize iteration counter
kiter = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Phase fraction of black phase
%f = 0.6;

% Pixel image dimensions
%N = 50;


% Target correlation function
% a0 = 10;         % conrols longer range correlations
% r0 = 5;         % controls the effective block size
% w = 2*pi/a0;
% s2_target = (1-f)^2 + f*(1-f)*exp(-r./r0).*sin(w*r)./(w*r);

% SEM image
I_targ = imageToMatrix('bin04.jpg');
I_targ = I_targ(801:900, 1401:1500);
[M,N] = size(I_targ);
f = sum(sum(I_targ))/(N^2);

% Max sampling length for calculation of correlation function
rsamp = N/2;
r = 0:rsamp;
r = r';

% Generate random initial image
I = generateImage(N,N,f);

% Save the initial image
I_init = I;

% Get the black and white positions from the initial image
[posB,posW] = getBWPositions(I);

% Calculate the target correlation function
s2_target = S2_pbc(I_targ,rsamp);

% Calculate the s2 correlation function and save the initial value
s2 = S2_pbc(I,rsamp);
s2_init = s2;
%s2 = fixS2(s2,s2_target,floor(0.8*rsamp));

% Initialize gross iteration counter at record vectors
C = 1;
moves = [];
Thist = [];
Ehist = [];

% Calculate the energy in the old configuration
E = calcEnergy(s2,s2_target);
Ehist(1) = E;
Thist(1) = T;
energies = zeros(L,1); 

% Simulated annealing loop
while E > Etol && T > Tmin
    
    % Loop through L times for a specific temperature
    for j = 1:L       
        
        % Swap pixels
        [Inew, posBnew, posWnew, r1, r2] = swapPixels(I,posB,posW);
        
        % Calculate change in s2 due to swap
        %ds2 = deltaS2_pbc(r1,r2,I,rsamp);
        %s2_new = s2 + ds2;
        % New correlation function
        s2_new = S2_pbc(Inew,rsamp);
        %s2_new = fixS2(s2_new,s2_target,floor(0.8*rsamp));
        
        % Calculate the new energy of the swapped image
        Enew = calcEnergy(s2_new, s2_target);
        
        % Calculate the change in energy for metropolis acceptance
        dE = Enew - E;
        
        % Metropolis algorithm
        [move, p, R] = myMetropolisMove(dE,T);

        % If we accept the new configuration, update the image, black/white
        % position vectors, and correlation function
        if  (move == 1)
            I = Inew;
            posB = posBnew;
            posW = posWnew;
            s2 = s2_new;
            E = Enew;
        end       
        
        % Update gross iteration counter
        C = C+1;
        moves(C) = move;
        energies(j) = E;
    end
    
% Update the iteration counter
kiter = kiter + 1;
% Record the temperature 
Thist(kiter) = T;
% Average E over the L trials for each temperature and record
Ehist(kiter) = sum(energies)/L;
% Reduce the temperature for next iteration
T = lambda*T;
end
%% Plot Correlation Functions
figure(1)
hold on
plot(r,s2,'bs')
plot(r,s2_target,'k-','Linewidth',2);
grid on
ll = legend('Reconstruction','Target');    set(ll,'Fontsize',12,'Location','Best');
xlabel('Distance (pixels)','Fontsize',12);
ylabel('Probability','Fontsize',12);
plot(r,s2_new,'b:');
title(sprintf('SEM Image Reconstruction (f = %4.2f)', f),'Fontsize',14);
hold off
suptitle('Correlation Functions')
%% Plot Optimization History
figure(2)
subplot(2,2,1)
plot(Thist)
xlabel('k')
ylabel('T')
subplot(2,2,2)
plot(moves,'bo')
xlabel('C')
ylabel('Acceptance')
subplot(2,2,[3,4])
plot(Ehist,'r-.')
xlabel('k')
ylabel('Energy')
suptitle('Optimization History')
%% Plot Image Matrices
figure(3)
subplot(1,2,1)
imagesc(I_init)
title('Initial')
subplot(1,2,2)
imagesc(I)
title('Final')
suptitle('Image Matrix Comparison')