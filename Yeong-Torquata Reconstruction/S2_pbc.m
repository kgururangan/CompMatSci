function [S2] = S2_pbc(I,rsamp)


% Get the image dimensions
[N,M] = size(I);

% select sampling range of line lengths to consider 
r = 0:rsamp;

% Initialize row and column 2-point correlation vectors
S2_row = zeros(rsamp+1,1);
S2_col = zeros(rsamp+1,1);

% Loop through the range of sample sizees
for jj = 1:length(r)
    
    % Loop through each row or column of the image
    for i = 1:N
        
        % Store row and column
        R = I(i,:);
        C = I(:,i);
        
        % Move a line of length r(jj) through the row or column by 1 unit
        for j = 1:N
            
            % Count the number of times a line of length r has both ends in
            % black phase 
            
            leftend = j;
            rightend = j+r(jj);
            
            % Periodic boundary conditions
            if (rightend > N)
                rightend = rightend - N;
            end
            
            if ( R(leftend) + R(rightend) == 2)
                    
                S2_row(jj) = S2_row(jj) + 1;
                
            end
            
            if (C(leftend) + C(rightend) == 2 )
                
                S2_col(jj) = S2_col(jj) + 1;
                
            end
            
            
        end
    end
    
end

% Total correlation function is average for isotropic medium
% Normalize over number of samples 
S2 = 0.5*(S2_row + S2_col)./N^2;
S2(end) = S2(1)^2;


end
