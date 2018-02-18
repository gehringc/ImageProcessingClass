clear all;
close all;

B = cell(1,7);
B{1} = imread('B1.bmp');
B{2} = imread('B2.bmp');
B{3} = imread('B3.bmp');
B{4} = imread('B4.bmp');
B{5} = imread('B5.bmp');
B{6} = imread('B6.bmp');
B{7} = imread('B7.bmp');

V = 255*[0.1 0 0.995; 
    0.8 0 0.6; 
    0.707 0 0.707; 
    -0.707 0 0.707;
    0 0.707 0.707; 
    0.56568 0.56568 0.6;
    -0.56568 0.56568 0.6;];
pinV = pinv(V);
I = [];
g = [];
G = [];
N = [];
albedo = [];
p = [];
q = [];
f = [];
width = size(B{1});
width = width(1);
for i=1:length(B{1})
    for j = 1:width
        for k = 1:7
            b = B{k};
            I(k,1) = b(j,i);            
        end
        G(j,i,:) = pinV*I;
        g = [G(j,i,1) G(j,i,2) G(j,i,3)];
        albedo(j,i) = norm(g);
        N(j,i,:) = g/albedo(j,i);
        p(j,i) = N(j,i,1)/N(j,i,3);
        q(j,i) = N(j,i,2)/N(j,i,3);
    end
end
f(1,1) = 0;
for i=2:340
    f(i,1) = f(i-1,1) + q(i,1);
end
for i=1:340
    for j=2:360
        f(i,j) = f(i,j-1) + p(i,j);
    end
end
figure(1);imshow(albedo);
figure(2);
Nx = N(:,:,1);
Ny = N(:,:,2);
Nz = N(:,:,3);
for i=1:5:340
    for j=1:5:360
        Nx2(i,j) = Nx(i,j);
        Ny2(i,j) = Ny(i,j);
        Nz2(i,j) = Nz(i,j);
        f2(i,j) = f(i,j);
    end
end
quiver3(f2,Nx2,Ny2,Nz2,5)
figure(3)
surf(f)
figure(4)
colormap(gray)
surf(f,albedo)