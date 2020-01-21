%% 5. Depth map computation with local methods
disp('... Runing exercice 5') 

% Data images: '0001_rectified_s.png','0002_rectified_s.png'

% Test the functions implemented in the previous section with the facade
% images. Try different matching costs and window sizes and comment the
% results.
% Notice that in this new data the minimum and maximum disparities may
% change.

imgLeft = imread('Data/0001_rectified_s.png');
imgRight = imread('Data/0002_rectified_s.png');
grayLeft = double(rgb2gray(imgLeft));
grayRight = double(rgb2gray(imgRight));
minDisp = 0;  
maxDisp = 16; 
gauss = 1;
matchingCost = 'SSD';  

w = 3;     
disparity = stereo_computation(grayLeft,grayRight,minDisp,mxd,w,matchingCost);
figure; imshow(uint8(disparity)*16);title('windowsize=3x3'); 

w = 9;    
disparity = stereo_computation(grayLeft,grayRight,minDisp,mxd,w,matchingCost);
figure; imshow(uint8(disparity)*16);title('windowsize=9x9');

w = 21;     
disparity = stereo_computation(grayLeft,grayRight,minDisp,mxd,w,matchingCost);
figure; imshow(uint8(disparity)*16); title('windowsize=21x21');

w = 31;     
disparity = stereo_computation(grayLeft,grayRight,minDisp,mxd,w,matchingCost);
figure; imshow(uint8(disparity)*16); title('windowsize=31x31');


matchingCost = 'NCC';  

w = 3;     
disparity = stereo_computation(grayLeft,grayRight,minDisp,mxd,w,matchingCost);
figure; imshow(uint8(disparity)*16);title('windowsize=3x3'); 

w = 9;    
disparity = stereo_computation(grayLeft,grayRight,minDisp,mxd,w,matchingCost);
figure; imshow(uint8(disparity)*16);title('windowsize=9x9');

w = 21;     
disparity = stereo_computation(grayLeft,grayRight,minDisp,mxd,w,matchingCost);
figure; imshow(uint8(disparity)*16); title('windowsize=21x21');

w = 31;     
disparity = stereo_computation(grayLeft,grayRight,minDisp,mxd,w,matchingCost);
figure; imshow(uint8(disparity)*16); title('windowsize=31x31');


