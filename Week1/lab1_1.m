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
I=imread('Data/0005_s.png'); % we have to be in the proper folder

% ToDo: generate a matrix H which produces a similarity transformation
% H is 3x3 metirx  H = [sR t; 0 1]
% R is a rotation matrix (orthogonal matrix) = [ cos(theta) - sin(theta);
% cos(theta) sin(theta)]
% t is a translation vector
t = [0;0];
s = 0.5; %isotropic scaling factor for scalated
theta = 15; %the orientation angle for rotated
H = [(s*cosd(theta)) (-s*sind(theta)) t_x;
     (s*sind(theta)) (s*cosd(theta)) t_y;
     0 0 1];

I2 = apply_H(I, H);
figure; imshow(I);title('original');
figure; imshow(uint8(I2));title('Similarities');


%% 1.2. Affinities

% ToDo: generate a matrix H which produces an affine transformation
% H is 3x3 H = [A t;0 1]
% A is a non-singular 2 × 2 matrix let A = [0.5 1; 1 0]
% t is a translation vector
s = 0.5;
theta = pi/4;
t = [0;0];
A = [0 1;1 0];
H = [0 1 t(1);
     1 0 t(2);
     0 0 1];

I2 = apply_H(I, H);
figure; imshow(uint8(I2));title(' Affinities');

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation
% SVD U-left singular vetors,V-right singular vetors D-diagonal metrix
[U,D,V] = svd(A);
Rtheta = U*V';
Rphi = V';
rotation_1 = [Rtheta(1,:) 0;Rtheta(2,:) 0;0 0 1];
rotation_2 = [Rphi(1,:) 0;Rphi(2,:) 0;0 0 1];
% scale
scale = [D(1,1) 0 0;0 D(2,2) 0;0 0 1];
% translation
translation = [1 0 t(1);0 1 t(2);0 0 1];

H_decompose= translation*(rotation_1*rotation_2*scale*rotation_2');

% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above
if(sum(round(H_decompose(:)-H(:)))) == 0
    disp('matrix H & H_decompose are equal');
else
    disp('matrix H & H_decompose are not equal');
end
% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before
I2_decompose = apply_H(I, H_decompose);
figure; imshow(uint8(I2_decompose));title('affine_decompose');


%% 1.3 Projective transformations (homographies)

% ToDo: generate a matrix H which produces a projective transformation
% The projective transformation Hp can be decomposed as Hs*Ha*Hp
t = [0;0];
s = 0.5; %isotropic scaling factor for scalated
theta = 15; %the orientation angle for rotated
% Hs similarity transformation
Hs = [(s*cosd(theta)) (-s*sind(theta)) t(1);
     (s*sind(theta)) (s*cosd(theta)) t(2);
     0 0 1];
% Ha affinity transformation 
A = [0 1;1 0];
Ha = [0 1 t(1);
     1 0 t(2);
     0 0 1];
% Hp projective transformation
% Hp=[I 0;vt v]
I = [1 0;0 1];
v = [0.00005,0.0005];
Hp = [I(1,:) t(1);I(2,:) t(2);v 0.7];

H = Hs*Ha*Hp;

I2 = apply_H(I, H);
figure; imshow(uint8(I2));title('projective');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification


% choose the image points
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
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

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points


% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% ToDo: compute the homography that affinely rectifies the image

I2 = apply_H(I, H);
figure; imshow(uint8(I2));

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4

% show the transformed lines in the transformed image
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that are used in every step.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)



