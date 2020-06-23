%% Simulated Annealing Reconstruction (Energy Threshold Scheme)

clear all
clc 
close all

% Phase fraction of black phase
f = 0.5;

% Pixel image dimensions
N = 50;

% Max sampling length for calculation of correlation function
rsamp = N/2;
r = 0:rsamp;
r = r';

% Initialize threshold energy
Eth = 1e6;

% Annealing schedule parameter
lambda = 0.99;

% Number of trials we want to take for each temperature
L = 10000;

% Threshold energy difference between target and image
Eth_min = 1e-9;
%dE = 1;
kiter = 0;

% Target correlation function
a0 = 5;
r0 = 4;
w = 2*pi/a0;
s2_target = (1-f)^2 + f*(1-f)*exp(-r./r0).*sin(w*r)./(w*r);
s2_target(1) = f;

% Generate image
I = generateImage(N,N,f);

% Save the initial image
I_init = I;

% Get the black and white positions from the initial image
[posB,posW] = getBWPositions(I);

% Calculate the s2 correlation function and save the initial value
s2 = S2_pbc(I,rsamp);
s2_init = s2;

% Initialize gross iteration counter and record vectors
C = 1;
thresh = [];
thresh(1) = Eth;
moves = [];
energies = [];

% Simulated annealing loop
while Eth > Eth_min 
    
    % Loop through L times for a specific temperature
    for j = 1:L
        
        % Calculate the energy in the old configuration
        Eold = calcEnergy(s2,s2_target);
        
        % Swap pixels
        [Inew, posBnew, posWnew, r1, r2] = swapPixels(I,posB,posW);
        
        % Calculate change in s2 due to swap
        ds2 = deltaS2_pbc(r1,r2,I,rsamp);
        
        % Calculate the new energy of the swapped image
        Enew = calcEnergy(s2 + ds2, s2_target);
        
        % Calculate the change in energy for metropolis acceptance
        dE = Enew - Eold;
        
        % Metropolis algorithm
        move = myEthAcceptance(dE,Eth);
        % Record move
        moves(C) = move;
        % Record energy
        energies(C) = Eold;
        
        % If we accept the new configuration, update the image and black/white
        % position vectors
        if  (move == 1)
            I = Inew;
            posB = posBnew;
            posW = posWnew;
            % If metropolis accepted, replace old energy with new
            energies(C) = Enew;
        end
        
        % Update gross iteration counter
        C = C+1;
        
    end

% Update iteration counter
kiter = kiter + 1;
% Record threshold energy
thresh(kiter) = Eth;
% Reduce the threshold energy
Eth = lambda*Eth;
end

%% Plot correlation functions
s2_final = S2_pbc(I,rsamp);
figure(1)
hold on
plot(r,s2_final,'ro',r,s2_target,'go',r,s2_init,'bo');
grid on
ll = legend('Final','Target','Initial');    set(ll,'Fontsize',12,'Location','Best');
xlabel('Distance (pixels)','Fontsize',12);
ylabel('Probability','Fontsize',12);
title(sprintf('k = %4.2f, f = %4.2f',kiter, f));
plot(r,s2_final,'r-',r,s2_target,'g-',r,s2_init,'b-');
hold off
%% Plot opthist
figure(2)
subplot(2,2,1)
plot(thresh)
xlabel('k')
ylabel('E_{th}')
subplot(2,2,2)
plot(moves,'b*')
xlabel('C')
ylabel('Acceptance')
subplot(2,2,[3,4])
plot(energies,'r*')
xlabel('C')
ylabel('E')
%% Plot image matrices
figure(3)
subplot(1,2,1)
imagesc(I_init)
subplot(1,2,2)
imagesc(I)
