function [dL] = deltaL_pbc(n,m,I)

% Get the image dimensions
[N,M] = size(I);

% Convert the positions of the switched pixels from linear indices to
% subscripts
[row1, col1] = ind2sub([N,M],n);
[row2, col2] = ind2sub([N,M],m);

% Extract relevant row and column vectors
R1 = I(row1,:);
R2 = I(row2,:);
C1 = I(:,col1);
C2 = I(:,col2);

% select sampling range for all chord length to consider 
rsamp = N;
r = 0:rsamp;

% Initialize row and column 2-point lineal path correlation vectors
dL_row = zeros(rsamp+1,1);
dL_col = zeros(rsamp+1,1);

% Sample over all chord lengths
for jj = 1:length(r)
    
        % Flag vectors store the positions of where the row has a 1
        flag_row1 = zeros(1,N);
        flag_col1 = flag_row1;
        flag_row2 = flag_row1;
        flag_col2 = flag_row1;
        
        % Loop through the entries of each row and column and record the
        % positions of 1 in flag
        for j = 1:N
                                         
            if (R1(j) == 1)
                flag_row1(j) = j;
            end
        
            if (C1(j) == 1)
                flag_col1(j) = j;
            end
            
            if (R2(j) == 1)
                flag_row2(j) = j;
            end
            
            if (C2(j) == 1)
                flag_col2(j) = j;
            end
        
        end
        
            % Calculate chord lengths held in flag
            % +1 is to account for matlab array indexing
            chord_r1 = abs(diff(flag_row1));
            chord_r2 = abs(diff(flag_row2));
            chord_c1 = abs(diff(flag_col1));
            chord_c2 = abs(diff(flag_col2));
        
        % Look through the chord lengths
        for k = 1:N-1           
            
            
            % Update L_row or Col by the formula 
            % L(z) = (chordlength - z)/N if 0 <= z <= chordlength
            %         0, otherwise
            
            % Look through the two rows
            
           if chord_r1(k) >= r(jj)
               
                dL_row(jj) = dL_row(jj)  +  (chord_r1(k) - r(jj))/N;
                
                % periodic bc
                if chord_r1(k) > flag_row1(1) + N - flag_row1(end)
                    % x2 is to account for the chord from last -> first and
                    % first -> last; THIS IS IMPORTANT!
                    dL_row(jj) = dL_row(jj) + 2*(chord_r1(k) - r(jj))/N;
                end
                    
                
           end
            
           if chord_r2(k) >= r(jj)
               
                dL_row(jj) = dL_row(jj)  +  (chord_r2(k) - r(jj))/N;
                 
                % periodic BC
                if chord_r2(k) > flag_row2(1) + N - flag_row2(end)
                    % x2 is to account for the chord from last -> first and
                    % first -> last; THIS IS IMPORTANT!
                    dL_row(jj) = dL_row(jj) + 2*(chord_r2(k) - r(jj))/N;
                end
                    
                
            end
            
            % The same for the two columns
            
            if chord_c1(k) >= r(jj)
                dL_col(jj) = dL_col(jj) + (chord_c1(k) - r(jj))/N;
                
                % Periodic b.c. x2 for degeneracy of last chord
                if chord_c1(k) > flag_col1(1) + N - flag_col1(end)
                    dL_col(jj) = dL_col(jj) + 2*(chord_c1(k) - r(jj))/N;
                end
            end
            

            if chord_c2(k) >= r(jj)
                dL_col(jj) = dL_col(jj) + (chord_c2(k) - r(jj))/N;
                
                % Periodic b.c. x2 for degeneracy of last chord
                if chord_c2(k) > flag_col2(1) + N - flag_col2(end)
                    dL_col(jj) = dL_col(jj) + 2*(chord_c2(k) - r(jj))/N;
                end
            end
            
            
        end
        
       
 

end

% Average and return lineal path correlation function
dL = 0.5*(dL_row + dL_col)./N^2;

end