clear all
close all
clc

f = 0.5;

N = 100;
b = 5;
T = 1000;
rsamp = N/2;

Ib = blockImage(f, [b,b], [N,N]);
I = generateImage(N,N,f);

s2 = S2_pbc(I,rsamp);

s2_target = S2_pbc(Ib,rsamp);

Eold = calcEnergy(s2,s2_target,rsamp);

% Get the black and white positions from the initial image
[posB,posW] = getBWPositions(I);

% Swap pixels 
[Inew, posBnew, posWnew, r1, r2] = swapPixels(I,posW,posW);
% Calculate change in s2 due to swap
ds2 = deltaS2_pbc(r1,r2,I,rsamp);
% Calculate the new energy of the swapped image
Enew = calcEnergy(s2 + ds2, s2_target, rsamp);
% Calculate the change in energy for metropolis acceptance
dE = Enew - Eold;
% Metropolis algorithm 
[move, p, R] = myMetropolisMove(dE,T);