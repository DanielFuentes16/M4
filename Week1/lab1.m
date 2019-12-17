%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 1: Image rectification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying image transformations

% ToDo: create the function  "apply_H" that gets as input a homography and
% an image and returns the image transformed by the homography.
% The size of the transformed image has to be automatically set so as to 
% contain the whole transformed image.
% At some point you will need to interpolate the image values at some points,
% you may use the Matlab function "interp2" for that.


%% 1.1. Similarities
close all;
I=imread('Data/0005_s.png'); % we have to be in the proper folder
figure; imshow(I);


% ToDo: generate a matrix H which produces a similarity transformation
angle = pi/4;
R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
t = [10; 0];
s = 2;

H = zeros(3,3);
H(1:2,1:2) = s*R;
H(1:2,3) = t;
H(3,3) = 1;
[I2, x_min, y_min] = apply_H(I, H);
figure; imshow(uint8(I2));



%% 1.2. Affinities

 %ToDo: generate a matrix H which produces an affine transformation
 close all;
 A = [1 0.5; 0.2 1]
 t = [10; 0]
 
 H = zeros(3,3);
 H(1:2,1:2) = A;
 H(1:2,3) = t; 
 H(3,3) = 1
 [I2, x_min, y_min] = apply_H(I, H);
  figure; imshow(uint8(I2));
 
% % ToDo: decompose the affinity in four transformations: two
% % rotations, a scale, and a translation
% 
[U,S,V] = svd(H(1:2,1:2));

s1 = S(1,1);
s2 = S(2,2);

R_phi = V';
phi_c = acos(R_phi(1,1));
phi_s = asin(R_phi(2,1));
if phi_s >= 0 % first or second quadrant
   phi = phi_c
elseif phi_c <= pi/2 % first or 4th quadrant
    phi = phi_s
else %3rd quadrant
    phi = 2*pi - phi_c
end

R_theta = U*V';
theta_c = acos(R_theta(1,1));
theta_s = asin(R_theta(2,1));
if theta_s >= 0 % first or second quadrant
    theta = theta_c
elseif theta_c <= pi/2 % first or 4th quadrant
    theta = theta_s
else %3rd quadrant
    theta = 2*pi - theta_c
end

Rotation_phi = zeros(3,3);
Rotation_phi(1:2,1:2) = R_phi;
Rotation_phi(3,3) = 1

Rotation_theta = zeros(3,3);
Rotation_theta(1:2,1:2) = R_theta;
Rotation_theta(3,3) = 1

Scale = zeros(3,3);
Scale(1:2,1:2) = S;
Scale(3,3) = 1

Translation = zeros(3,3);
Translation(1,1) = 1;
Translation(2,2) = 1;
Translation(3,3) = 1;
Translation(1:2,3) = H(1:2,3)
% 
% 
% % ToDo: verify that the product of the four previous transformations
% % produces the same matrix H as above
% 
H_recomposed = Translation*Rotation_theta*inv(Rotation_phi)*Scale*Rotation_phi
% 
% % ToDo: verify that the proper sequence of the four previous
% % transformations over the image I produces the same image I2 as before
% 
[I3, x_min, y_min] = apply_H(I, Rotation_phi);
figure; imshow(uint8(I3));

[I4, x_min, y_min] = apply_H(I3, Scale);
figure; imshow(uint8(I4));

[I5, x_min, y_min] = apply_H(I4, inv(Rotation_phi));
figure; imshow(uint8(I5)); 

[I6, x_min, y_min] = apply_H(I5, Rotation_theta);
figure; imshow(uint8(I6)); 

[I7, x_min, y_min] = apply_H(I6, Translation);
figure; imshow(uint8(I7)); 
% 
% 
% 
%% 1.3 Projective transformations (homographies)

% ToDo: generate a matrix H which produces a projective transformation
H(3,1:2) = [0 0.001];
I2 = apply_H(I, H);
figure; imshow(uint8(I2));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification
% 
% 
% % choose the image points
 I=imread('Data/0000_s.png');
 A = load('Data/0000_s_info_lines.txt');
% 
% % indices of lines
 i = 424;
 p1 = [A(i,1) A(i,2) 1]';
 p2 = [A(i,3) A(i,4) 1]';
 i = 240;
 p3 = [A(i,1) A(i,2) 1]';
 p4 = [A(i,3) A(i,4) 1]';
 i = 712;
 p5 = [A(i,1) A(i,2) 1]';
 p6 = [A(i,3) A(i,4) 1]';
 i = 565;
 p7 = [A(i,1) A(i,2) 1]';
 p8 = [A(i,3) A(i,4) 1]';
% 
% % ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
 [l1, l2, l3, l4] = compute_lines(p1, p2, p3, p4, p5, p6, p7, p8);
% 
% % show the chosen lines in the image
faar = figure();
figure(faar);
imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'b');
% % ToDo: compute the homography that affinely rectifies the image
 H = get_affinity(l1, l2, l3, l4);
 I2 = apply_H(I, H);
 
% % 
% % ToDo: compute the transformed lines lr1, lr2, lr3, lr4
 H_inv_t = inv(H');
 lr1 = H_inv_t*l1';
 lr2 = H_inv_t*l2';
 lr3 = H_inv_t*l3';
 lr4 = H_inv_t*l4';
 
% % 
% % show the transformed lines in the transformed image
fa = figure();
figure(fa)
imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2) , 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2) , 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2) , 'g');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2) , 'g');


% % show the transformed lines in the transformed image
fb = figure();
figure(fb)
imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2) , 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2) , 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2) , 'g');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2) , 'g');

% 
% 
% 
% % % ToDo: to evaluate the results, compute the angle between the different pair 
% % % of lines before and after the image transformation
% 
compute_angle(l1,l2)
compute_angle(l3,l4)

compute_angle(lr1,lr2)
compute_angle(lr3,lr4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 3. Metric Rectification
% 
% %% 3.1 Metric rectification after the affine rectification (stratified solution)
% 
% % ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
% %       As evaluation method you can display the images (before and after
% %       the metric rectification) with the chosen lines printed on it.
% %       Compute also the angles between the pair of lines before and after
% %       rectification.
% 
lmr1 = lr1;
lmr2 = lr3;
lmr3 = cross(cross(lr1,lr3),cross(lr2,lr4));
lmr4 = cross(cross(lr2,lr3),cross(lr1,lr4));
hold off;
fbm = figure();
figure(fbm)
imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lmr1(1)*t + lmr1(3)) / lmr1(2), 'y');
plot(t, -(lmr2(1)*t + lmr2(3)) / lmr2(2), 'y');
plot(t, -(lmr3(1)*t + lmr3(3)) / lmr3(2), 'c');
plot(t, -(lmr4(1)*t + lmr4(3)) / lmr4(2), 'c');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'k');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'k');
% 
% %1,2. get (s1,s2,s3) for the system of equations defined by the two
% %orthogonal pairs
syms s1 s2 s3
eqn1 = lmr1(1)*lmr2(1)*s1 + (lmr1(1)*lmr2(2)+lmr1(2)*lmr2(1))*s2 + lmr1(2)*lmr2(2)*s3 == 0;
eqn2 = lmr3(1)*lmr4(1)*s1 + (lmr3(1)*lmr4(2)+lmr3(2)*lmr4(1))*s2 + lmr3(2)*lmr4(2)*s3 == 0;
eqn3 = s3 == 1;
X = solve(eqn1, eqn2, eqn3);
S = [X.s1 X.s2; X.s2 X.s3]
% 
% %3. Knowing S, use the Cholesky decomposition of S to compute an upper
% %triangular matrix K such that S = KK T
K = chol(S);
% 
% %4,5. Calculate Ha<-s
Ha_s = zeros(3, 3);
Ha_s(1:2,1:2) = K;
Ha_s(3,3) = 1;

I_metric_rectified = apply_H(I2, Ha_s);
% 
% % compute the transformed lines lamr1, lamr2, lamr3, lamr4
 H_mr = inv(Ha_s');
 lamr1 = H_mr*lmr1;
 lamr2 = H_mr*lmr2;
 lamr3 = H_mr*lmr3;
 lamr4 = H_mr*lmr4;
% % show the transformed lines in the transformed image
famr = figure();
figure(famr)
imshow(uint8(I_metric_rectified));
hold on;
t=1:0.1:1000;
plot(t, -(lamr1(1)*t + lamr1(3)) / lamr1(2), 'y');
plot(t, -(lamr2(1)*t + lamr2(3)) / lamr2(2), 'y');
plot(t, -(lamr3(1)*t + lamr3(3)) / lamr3(2), 'c');
plot(t, -(lamr4(1)*t + lamr4(3)) / lamr4(2), 'c');
% 
% %       Compute also the angles between the pair of lines before and after
% %       rectification.
% 
pi/2
compute_angle(lamr1, lamr2)
compute_angle(lamr3, lamr4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that are used in every step.
clear;
I = imread('Data/0001_s.png');
A = load('Data/0001_s_info_lines.txt');

% Crop the initial image so that only the left facade
I = I(:,1:470,:); 

% indices of lines
i = 614; 
i = 188; 
p1 = [A(i,1) A(i,2) 1]'; 
p2 = [A(i,3) A(i,4) 1]'; 
i = 159; 
p3 = [A(i,1) A(i,2) 1]'; 
p4 = [A(i,3) A(i,4) 1]'; 
i = 645; 
p5 = [A(i,1) A(i,2) 1]'; 
p6 = [A(i,3) A(i,4) 1]'; 
i = 541; 
p7 = [A(i,1) A(i,2) 1]'; 
p8 = [A(i,3) A(i,4) 1]'; 

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
    point1 = [p1(1) p1(2) 1];
    point2 = [p2(1) p2(2) 1];
    point3 = [p3(1) p3(2) 1];
    point4 = [p4(1) p4(2) 1];
    point5 = [p5(1) p5(2) 1];
    point6 = [p6(1) p6(2) 1];
    point7 = [p7(1) p7(2) 1];
    point8 = [p8(1) p8(2) 1];

    l1 = cross(point1, point2);
    l2 = cross(point3, point4);
    l3 = cross(point5, point6);
    l4 = cross(point7, point8);


% show lines 
figure(1); imshow(I); title('Original');
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'r');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'r');

% ToDo: compute the homography that affinely rectifies the image
% compute vanishing points where lines cross at ininity
v1 = cross(l1, l2);
v2 = cross(l3, l4);

% compute line that passes through vanishing points: line at infinite 
l_inf = cross(v1, v2); 
l_inf = l_inf / l_inf(3);

% define H for affine rectification and apply to image
H = [1 0 0; 0 1 0; l_inf(1) l_inf(2) 1];
I_trans = apply_H(I, H);
figure(2);imshow(uint8(I_trans)); title('Affine rectification');

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
% by using l' = H^(-T) * l
vanishing_point1 = cross(l1,l2);
vanishing_point1 = vanishing_point1/vanishing_point1(3);
vanishing_point2 = cross(l3,l4);
vanishing_point2 = vanishing_point2/ vanishing_point2(3);
vanishing_line = cross(vanishing_point1,vanishing_point2);

H = [1 0 0;
     0 1 0; 
     vanishing_line(1)/vanishing_line(3) vanishing_line(2)/vanishing_line(3) 1];
% transformed points
HP1 = H*p1; HP2 = H*p2; HP3 = H*p3; HP4 = H*p4;
HP5 = H*p5; HP6 = H*p6; HP7 = H*p7; HP8 = H*p8;
 
HP1 = HP1/HP1(3); HP2 = HP2/HP2(3); HP3 = HP3/HP3(3); HP4 = HP4/HP4(3);
HP5 = HP5/HP5(3); HP6 = HP6/HP6(3); HP7 = HP7/HP7(3); HP8 = HP8/HP8(3);

% transformed lines
lr1 = cross(HP1,HP2);
lr2 = cross(HP3,HP4);
lr3 = cross(HP5,HP6);
lr4 = cross(HP7,HP8);

% show the transformed lines in the transformed image
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'r');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'r');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation
% calculate slopes of original lines
s1 = l1(1) / l1(2); 
s2 = l2(1) / l2(2);
s3 = l3(1) / l3(2);
s4 = l4(1) / l4(2);
% calculate slope of rectified lines
sr1 = lr1(1) / lr1(2); 
sr2 = lr2(1) / lr2(2);
sr3 = lr3(1) / lr3(2);
sr4 = lr4(1) / lr4(2);

% angles between "parallel" original lines
a12 = rad2deg( atan(s1) - atan(s2) ); % a12 = 0.0992
a34 = rad2deg( atan(s3) - atan(s4) ); % a34 = -1.3435
sprintf('Angles between yellow and red lines %f, %f', a12, a34)
% angles of really parallel rectified lines
ar12 = rad2deg( atan(sr1) - atan(sr2) ); % ar12 = -7.9514e-16 (almost 0)
ar34 = rad2deg( atan(sr3) - atan(sr4) ); % ar34 = 0
sprintf('Angles between yellow and red lines (after rectification) %f, %f', ...
  ar12, ar34)

% Orthogonal pair of lines (l1,r1) and (l2,m2)
l1 = lr1;
m1 = lr3;
l2 = lr2;
m2 = lr4;

A = [ l1(1)*m1(1), l1(1)*m1(2)+l1(2)*m1(1), l1(2)*m1(2);
      l2(1)*m2(1), l2(1)*m2(2)+l2(2)*m2(1), l2(2)*m2(2) ];
 
% A = [l1(1)*m1(1),   l1(1)*m1(2)+l1(2)*m1(1),    l1(2)*m1(2)];

s = null(A);
S = [ s(1), s(2); s(2), s(3) ] / s(3);

% Apply Cholesky factorization
K = chol(S); 
H = eye(3);
K = inv(K); 
H(1:2,1:2) = K;
H = H';

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
% l' = H^(-T) * l
lr1 = H' \ l1; lr1 = lr1 / lr1(3);
lr2 = H' \ l2; lr2 = lr2 / lr2(3);
m1r = H' \ m1; m1r = m1r / m1r(3);
m2r = H' \ m2; m2r = m2r / m2r(3);

% show transformed lines
I3 = apply_H(uint8(I_trans), H);
figure(3);imshow(uint8(I3)); title('Metric rectification');
hold on;
t=1:0.1:10000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(m1r(1)*t + m1r(3)) / m1r(2), 'r');
plot(t, -(m2r(1)*t + m2r(3)) / m2r(2), 'r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)



