
function [dS2] = deltaS2_pbc(n,m,I,rsamp)

% Given a pixel swap at linear indices n and m, the change in the 2-point
% correlation function (ds2 = s2_new - s2_old) is given by the 2-point
% correlation function over the rows and columns of the swapped pixels

% get image dimensions
[N,M] = size(I);

% convert linear indices of swapped pixels to rows and cols
[row1,col1] = ind2sub([N,M],n);
[row2,col2] = ind2sub([N,M],m);

% specify sampling range
r = 0:rsamp;

% initialize row and column 2-point correlation vectors
dS2_row = zeros(rsamp+1,1);
dS2_col = zeros(rsamp+1,1);

% sample over all line lengths (increase the line length you are sampling
for jj = 1:length(r)
    
    % extract relevant rows and columns
    R1 = I(row1,:);
    R2 = I(row2,:);
    C1 = I(:,col1);
    C2 = I(:,col2);
    
    % move the line by one pixel through the row or column each iteration
    for j = 1:N
    
        % check if endpoints of line lie in the black phase for each row
        % and column
        
            leftend = j;
            rightend = r(jj) + j;
            
            % Periodic boundary conditions
             if rightend > N
                rightend = rightend - N;
             end
                    
             if (R1(leftend) + R1(rightend) == 2)        
                dS2_row(jj) = dS2_row(jj) + 1;        
             end
    
             if (R2(leftend) + R2(rightend) == 2)        
                dS2_row(jj) = dS2_row(jj) + 1;
             end
    
             if (C1(leftend) + C1(rightend) == 2)        
                dS2_col(jj) = dS2_col(jj) + 1;        
             end
    
             if (C2(leftend) + C2(rightend) == 2)        
                dS2_col(jj) = dS2_col(jj) + 1;
             end
    
    end
    
    
end 


    % calculate total change in s2 as the average of the row and column changes
    % normalize the result by the number of sample sizes (N)
    dS2 = 0.5*(dS2_row + dS2_col)./N^2;
    %dS2(end) = 0;
    dS2(1) = 0;


end

