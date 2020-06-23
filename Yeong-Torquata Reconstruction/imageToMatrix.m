function [ M ] = imageToMatrix( img )
% This function converts a RGB image file into B/W with matrix of 1/0
%   img is the image to convert in the format of 'img'
RGB = imread(img);
figure
imshow(RGB)

I = rgb2gray(RGB);
figure
imshow(I)

BW = imbinarize(I,'global');
figure
imshow(BW)

M = int8(BW);
end

