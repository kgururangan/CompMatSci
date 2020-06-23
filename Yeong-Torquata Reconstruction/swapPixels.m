function [ Inew, posBnew, posWnew, b, w ] = swapPixels( I, posB, posW)

    % Function randomly swaps two pixels of different color (black or
    % white)
    
    % Input: posB - array of linear index positions for 1's in I
    %        posW - array of linear index positions for 0's in I
    %        I - original image
    % Output: I - swapped image

    Inew = I;
    posBnew = posB;
    posWnew = posW;
    
    % Randomly select an element from both posB and posW
    r1 = randi(length(posB),[1,1]);
    r2 = randi(length(posW),[1,1]);
    
    % Linear indices of pixels to be swapped
    b = posB(r1);   % index of black pixel
    w = posW(r2);   % index of white pixel

    % Switch the black and white pixels in the image matrix
    temp = I(w);
    Inew(w) = I(b);
    Inew(b) = temp;

    % Update the position vectors
    posBnew(r1) = w;
    posWnew(r2) = b;




end
