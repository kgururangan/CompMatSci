function [ B, W ] = getBWPositions( I )

% Input: I - image matrix of zeros and ones
% Output: B - vector containing linear index positions of where 1's are
%         W - vector containing linear index positiosn of where 0's are

[N,M] = size(I);

cb = 1;
cw = 1;

B = [];
W = [];

for i = 1:N*M
    
    if (I(i) == 1)
        
        B(cb) = i;
        cb = cb + 1;
        
    else
        
        W(cw) = i;
        cw = cw + 1;
        
    end
   


end

end
