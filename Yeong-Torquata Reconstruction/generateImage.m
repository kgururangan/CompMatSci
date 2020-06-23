function [ I ] = generateImage( N, M, f )

    % Generates an NxM matrix of 1's and 0's where 1 = phase of interest, 0
    % = matrix phase
    % The phase fraction of 1 is f where 0 < f < 1
    
    nblack = ceil(f*N*M);

    I = zeros(N,M);
    
    % Randomly select nblack unique points from the the list of 1:N*M
    R = randperm(N*M,nblack);
    
    I(R) = 1;
  
    

end

