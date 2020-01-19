function [result] = stereo_computation(leftImage, rightImage, minDisp, maxDisp, windowsSize, cost)
%STEREO_COMPUTATION Summary of this function goes here
%   Detailed explanation goes here
% The input parameters are 5:
% - left image
% - right image
% - minimum disparity
% - maximum disparity
% - window size (e.g. a value of 3 indicates a 3x3 window)
% - matching cost (the user may able to choose between SSD and NCC costs)

[w h] = size(leftImage);
halfWindows = fix(windowsSize/2);
result = ones(size(leftImage,1), size(leftImage,2));
    for i = 1 + halfWindows:w -halfWindows
        for j = 1 + halfWindows:h-halfWindows - maxDisp
            
            target = rightImage(i - halfWindows : i + halfWindows, j - halfWindows : j + halfWindows);
            count = 1;
            bestResult = 0;
            bestCost = inf;
            c =0;
            for countDis = minDisp : maxDisp
                %if  j + countDis - halfWindows > 0 && j + countDis + halfWindows < h
                    source = leftImage(i - halfWindows : i + halfWindows, j+ countDis - halfWindows : j + countDis + halfWindows);
                    if(cost == 'SSD')
                        %Minimize sum of square differences
                        c = sum(sum((w(:).*((target(:) - source(:)).^2))));
                                            
                    end
                    if (cost == 'NCC')
                         mean1 = mean(target(:));
                         mean2 = mean(source(:));

                        c = sum(w*(target-mean1).*(source-mean2))/(sqrt(sum(sum(w*(target-mean1).^2)))*sqrt(sum(sum(w*(source-mean2).^2))));
                    end
                    if (cost == 'BIL')
                        T = 40;
                        gammaC = 5;
                        p = [halfWindows+1,halfWindows+1];
                        wLeft = zeros(size(target));
                        wRight = zeros(size(source));
                        deltaCLeft = abs(target(p(1),p(2))-target); 
                        deltaCRight = abs(source(p(1),p(2))-source); 

                        wLeft = exp(-deltaCLeft/gammaC).*exp(-zeros((halfWindows*2)+1)/halfWindows);
                        wRight = exp(-deltaCRight/gammaC).*exp(-zeros((halfWindows*2)+1)/halfWindows);

                        w=wLeft.*wRight;
                        c = sum(sum((w(:).*((source(:) - target(:)).^2))));
                    end
                    %result(i , j, count) = c - 1;
                    %count = 1 + count;
                %end
                if c < bestCost
                    bestCost = c;
                    bestResult = countDis;
                end
            end
            result(i, j)=bestResult;
        end
    end
end