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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)



