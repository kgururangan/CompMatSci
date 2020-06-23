function [ image ] = blockImage( f,blocksize,dimension )
% This function takes in a fraction of 1s over 0s, blocksize of 1s, and
% dimension of the final binary matrix. This function outputs a binary matrix
% Blocksize and dimension are 1x2 matrices, and f is a fraction.

M = dimension(1);
N = dimension(2);
m = blocksize(1);
n = blocksize(2);

image = ones(dimension);
nZeros = f*M*N;
nBlocks = ceil(nZeros/(m*n));
IR = [];
IR(1) = 1;
IC = [];
IC(1) = 1;

for i = 1:floor((M-m)/m)
    IR(i+1) = IR(i) + m;
end
for j = 1:floor((N-n)/n)
    IC(j+1) = IC(j) + n;
end

c = 1;
output = zeros(1,2);

while c <= nBlocks
    flags = 0;
    new = [randi(length(IR)),randi(length(IC))];
    temp = [output; new];
    for k = 1:c
        if new == output(k,:)
            flags = flags + 1;
        end
    end
    if flags == 0;
    output = temp;
    c = c+1;
    else
        continue
    end
            
end

output(1,:) =[];

for ii = 1:nBlocks
    for jj = 0:m-1
        for kk = 0:n-1
            image(IR(output(ii,1))+jj,IC(output(ii,2))+kk) = 0;
        end
    end

end

% Flip the colors
image = -(image - 1);

end
