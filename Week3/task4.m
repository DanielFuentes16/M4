%% 4. OPTIONAL:  with your own images

% 4.1 Take a set of images of a moving scene from different viewpoints at 
%     different time instants. At least two images have to be taken from
%     roughly the same location by the same camera.
%
% 4.2 Implement the first part (until line 16) of the Algorithm 1 of the 
%     Photo-sequencing paper with a selection of the detected dynamic
%     features. You may reuse the code generated for the previous question.
%

%     features. You may reuse the code generated for the previous question.

im1 = imread('Data/im1.jpg');
im2 = imread('Data/im2.jpg');
im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;
[points_1, desc_1] = sift(im1, 'Threshold', 0.01);
[points_2, desc_2] = sift(im2, 'Threshold', 0.01);
matches = siftmatch(desc_1, desc_2);
    
% p1 and p2 contain the homogeneous coordinates of the matches
p1 = [points_1(1:2, matches(1,:)); ones(1, length(matches))];
p2 = [points_2(1:2, matches(2,:)); ones(1, length(matches))];
    
figure;
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches, 'Stacking', 'v');
   
[fD1, fS1, fDk, fSk] = Classify_Dyn_Stat(p1, p2);                        
for k = N-2                                                               
    [matches, p1, pk] = matchOpcional(imread('Data/im1.jpg'), imread(strcat('Data/im',int2str(K+2),'.jpg'))); 
    [fDk, fSk] = Classify_Dyn_Stat2(p1, pk);                    
    [F, inliers] = ransac_fundamental_matrix(fDk, fSk, 2.0);                 
end                                                                        

for df =...;                                                               
    l== p1 * p2;                                                             
    for pk= ...;                                                            
        lk== ...*;                                                         
        pk == l*lk;                                                          
        a = ComputeAlpha(p1,p2, pk);                                       
    end                                                                     
    o = sort(a);                                                            
end                                                                         
