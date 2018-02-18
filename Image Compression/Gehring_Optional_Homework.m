clear all;
% RGB = imread('4900.jpg');
% R = rgb2gray(RGB);
% J = dct2(R);
% imshow(J), colormap(jet(64))
% J = J(2340/4:2340/1.5,4160/4:4160/1.5);
% J(abs(J) < 10) = 0;
% K = idct2(J);
% figure, imshow(K)

RGB = imread('4900.jpg');
I = rgb2gray(RGB);
imwrite(I,'4901.jpg');
I = double(rgb2gray(RGB));
dct = @(block_struct) dct2(block_struct.data);
idct = @(block_struct) idct2(block_struct.data);
B = blockproc(I,[8,8],dct);
mask = find(abs(B) < 1000);
B(mask) = zeros(size(mask));
K = blockproc(B, [8 8], idct)/255;
imshow(K);
imwrite(K,'4900_compressed.jpg');