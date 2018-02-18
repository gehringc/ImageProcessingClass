clear all; close all;
%This program takes points from a 2D image of a calibration grid and
%uses its coordinates on the grid to redraw the calibration grid in color
%% Load image and hand pick points
C = imread('project1_image1.jpg');
info = imfinfo('project1_image1.jpg');
% Points on grid
P = [5 0 5 1; 0 0 0 1; 2 2 0 1; 0 10 0 1; 0 6 6 1; 0 3 9 1; 8 0 3 1; 3 8 0 1];
% The points on the image that translate to the grid points (pixels)
p = [1907 1221; 1345 1804; 1495 1889; 690 2330; 961 1180; 1128 709; 2288 1477; 1331 2219];
zero = [0;0;0;0];
l = 0;
k = 1;
%% Create perspective projection matrix
% U is a matrix containing 2 linear equations for each of the 3D points and 
% the 2D vertices x, and y. So 16 equations
for i=1:length(P)
    for j=1:4
        U(k,j) = P(i,j);
    end
    for j=5:8
        U(k,j) = 0;
    end
    for j=1:4
        U(k,j+8) = -p(i,1)*P(i,j);
    end
    l = l + 4;
    k = k + 1;
    for j=1:4
        U(k,j) = 0;
    end
    for j=1:4
        U(k,j+4) = (P(i,j));
    end
    for j=1:4
        U(k,j+8) = -p(i,2)*P(i,j);
    end
    k = k + 1;
end
Ut = transpose(U);
%Use linear least squares to minimize solution
lambda = Ut*U;
%Find eigenvalues and eigenvectors to find minimizing solution
[V,D] = eig(lambda);
solution = D(1,1)*V(:,1);
m1 = [solution(1) solution(2) solution(3) solution(4)];
m2 = [solution(5) solution(6) solution(7) solution(8)];
m3 = [solution(9) solution(10) solution(11) solution(12)];
%Create 3x4 perspective projection matrix using solution
M = [m1;m2;m3];

%% Find extrinsic and intrinsic parameters

A = [M(:,1) M(:,2) M(:,3)];
b = M(:,4);
a1t = A(:,1);
a2t = A(:,2);
a3t = A(:,3);
epsilon = 1;
rho = epsilon/sqrt(a3t(1)^2 + a3t(2)^2 + a3t(3)^2);
r3 = rho*transpose(a3t);
x0 = rho^2*(transpose(a1t)*a3t);
y0 = rho^2*(transpose(a2t)*a3t);
a1xa3 = cross(transpose(a1t),transpose(a3t));
a2xa3 = cross(transpose(a2t),transpose(a3t));
cosTheta = -(a1xa3*transpose(a2xa3))/(sqrt(a1xa3(1)^2 + a1xa3(2)^2 + a1xa3(3)^2)*sqrt(a2xa3(1)^2 + a2xa3(2)^2 + a2xa3(3)^2));
theta = acos(cosTheta);
alpha = rho^2*sqrt(a1xa3(1)^2 + a1xa3(2)^2 + a1xa3(3)^2)*sin(theta);
beta = rho^2*sqrt(a2xa3(1)^2 + a2xa3(2)^2 + a2xa3(3)^2)*sin(theta);
r1 = (1/sqrt(a2xa3(1)^2 + a2xa3(2)^2 + a2xa3(3)^2))*a2xa3;
r2 = cross(r3,r1);
K = [alpha -alpha*cot(theta) x0; 0 beta/sin(theta) y0; 0 0 1];
t = K\b;
t = rho*t;

%% Using perspective projection matrix and parameters, replot grid in color

i = 1;
%plot the image
figure,imshow(C)
    hold on;
%plot each line of grid using plot function to plot a line between
%two points. 3D points are used and then projected using the equation
% p = (1/Z)*M*P
for q=1:9
    top(q,:) = [q 0 10 1];
    bottom(q,:) = [q 0 0 1];
    topp(q,:) = (1/top(1,3))*M*transpose(top(q,:));
    botp(q,:) = (1/0.0000000000001)*M*transpose(bottom(q,:));
    topp1(q,:) = [int64(topp(q,1)/topp(q,3)) int64(topp(q,2)/topp(q,3))];
    botp1(q,:) = [int64(botp(q,1)/botp(q,3)) int64(botp(q,2)/botp(q,3))];
    plot([botp1(q,1),topp1(q,1)],[botp1(q,2),topp1(q,2)],'Color','r','LineWidth',2)
end
for q=1:9
    top2(q,:) = [0 0 q 1];
    bottom2(q,:) = [10 0 q 1];
    topp2(q,:) = (1/top2(q,3))*M*transpose(top2(q,:));
    botp2(q,:) = (1/bottom2(q,3))*M*transpose(bottom2(q,:));
    topp12(q,:) = [int64(topp2(q,1)/topp2(q,3)) int64(topp2(q,2)/topp2(q,3))];
    botp12(q,:) = [int64(botp2(q,1)/botp2(q,3)) int64(botp2(q,2)/botp2(q,3))];
    plot([botp12(q,1),topp12(q,1)],[botp12(q,2),topp12(q,2)],'Color','r','LineWidth',2)
end
for q=1:9
    top3(q,:) = [0 q 10 1];
    bottom3(q,:) = [0 q 0 1];
    topp3(q,:) = (1/top3(1,3))*M*transpose(top3(q,:));
    botp3(q,:) = (1/0.0000000000001)*M*transpose(bottom3(q,:));
    topp13(q,:) = [int64(topp3(q,1)/topp3(q,3)) int64(topp3(q,2)/topp3(q,3))];
    botp13(q,:) = [int64(botp3(q,1)/botp3(q,3)) int64(botp3(q,2)/botp3(q,3))];
    plot([botp13(q,1),topp13(q,1)],[botp13(q,2),topp13(q,2)],'Color','b','LineWidth',2)
end
for q=1:9
    top4(q,:) = [0 0 q 1];
    bottom4(q,:) = [0 10 q 1];
    topp4(q,:) = (1/top4(q,3))*M*transpose(top4(q,:));
    botp4(q,:) = (1/bottom4(q,3))*M*transpose(bottom4(q,:));
    topp14(q,:) = [int64(topp4(q,1)/topp4(q,3)) int64(topp4(q,2)/topp4(q,3))];
    botp14(q,:) = [int64(botp4(q,1)/botp4(q,3)) int64(botp4(q,2)/botp4(q,3))];
    plot([botp14(q,1),topp14(q,1)],[botp14(q,2),topp14(q,2)],'Color','b','LineWidth',2)
end
for q=1:9
    top5(q,:) = [q 0 0 1];
    bottom5(q,:) = [q 10 0 1];
    topp5(q,:) = (1/0.0000000000001)*M*transpose(top5(q,:));
    botp5(q,:) = (1/0.0000000000001)*M*transpose(bottom5(q,:));
    topp15(q,:) = [int64(topp5(q,1)/topp5(q,3)) int64(topp5(q,2)/topp5(q,3))];
    botp15(q,:) = [int64(botp5(q,1)/botp5(q,3)) int64(botp5(q,2)/botp5(q,3))];
    plot([botp15(q,1),topp15(q,1)],[botp15(q,2),topp15(q,2)],'Color','g','LineWidth',2)
end
for q=1:9
    top6(q,:) = [0 q 0 1];
    bottom6(q,:) = [10 q 0 1];
    topp6(q,:) = (1/0.0000000000001)*M*transpose(top6(q,:));
    botp6(q,:) = (1/0.0000000000001)*M*transpose(bottom6(q,:));
    topp16(q,:) = [int64(topp6(q,1)/topp6(q,3)) int64(topp6(q,2)/topp6(q,3))];
    botp16(q,:) = [int64(botp6(q,1)/botp6(q,3)) int64(botp6(q,2)/botp6(q,3))];
    plot([botp16(q,1),topp16(q,1)],[botp16(q,2),topp16(q,2)],'Color','g','LineWidth',2)
end
