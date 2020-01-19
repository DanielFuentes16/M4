function disparity = stereo_computation(leftImage, rightImage, minDisp, maxDisp, windowsSize, cost, gauss)
%STEREO_COMPUTATION Summary of this function goes here
%   Detailed explanation goes here
% The input parameters are 5:
% - left image
% - right image
% - minimum disparity
% - maximum disparity
% - window size (e.g. a value of 3 indicates a 3x3 window)
% - matching cost (the user may able to choose between SSD and NCC costs)

    [rows, cols] = size(leftImage);
    disparity = zeros(size(leftImage));

    %Window distance from center
    windowDist = ceil(windowsSize/2)-1;
    
    deltaG = zeros((windowDist*2)+1);
    for i=-windowDist:windowDist
        for j=-windowDist:windowDist    
            deltaG(i+windowDist+1,j+windowDist+1) = sqrt((i^2)+(j^2));
        end
    end
        
    for i = windowDist + 1:rows - windowDist
        for j = windowDist + 1:cols - windowDist
            %Left region
            regionLeft = rightImage((i - windowDist):(i + windowDist),(j - windowDist):(j + windowDist));            
            
            %Set weigths, gauss = 1 for section 5
            w = 1/(prod(size(regionLeft)))*ones(size(regionLeft));
            %Optional sigma values  = 0.5, 0.7, 0.9
            if (gauss == 1)
                sigma = 0.7;
                w = fspecial('gaussian', size(regionLeft), sigma); 
            end    

            %Find initial and final position for disparity search 
            ini_idx = j + minDisp; 
            fin_idx = j + maxDisp;
            
            %Control that disparity search is not out of image range
            ini_s = max(ini_idx, 1 + windowDist);
            fin_s = min(cols-windowDist, fin_idx);
            
            %max/min values for matching costs
            maxValue = 0;
            minValue = Inf;
            
            %bestDisplay is disparity value with which the minimum (SSD) 
            %or maximum (NCC) cost value is found
            bestDisplay = 0;
            for  s = ini_s:fin_s
                if(abs(s-j) >= minDisp)
                    regionRight = leftImage((i - windowDist):(i + windowDist),(s - windowDist):(s + windowDist));
                    switch cost
                        
                        case 'SSD'
                            %Minimize sum of square differences
                            ssd = sum(sum((w(:).*((regionLeft(:) - regionRight(:)).^2))));
                            if ssd < minValue
                                minValue = ssd;
                                bestDisplay = s - ini_s;
                            end
                            
                        case 'NCC'
                            %Normalized Cross Correlation
                            mean1 = mean(regionLeft(:));
                            mean2 = mean(regionRight(:));

                            ncc = sum(w*(regionLeft-mean1).*(regionRight-mean2))/(sqrt(sum(sum(w*(regionLeft-mean1).^2)))*sqrt(sum(sum(w*(regionRight-mean2).^2))));
                            if ncc > maxValue
                                maxValue = ncc;
                                bestDisplay =  s - ini_s;
                            end
                            
                        case 'bilateral'
                            T = 40;
                            gammaC = 5;
                            p = [windowDist+1,windowDist+1];
                            wLeft = zeros(size(regionLeft));
                            wRight = zeros(size(regionRight));
                            deltaCLeft = abs(regionLeft(p(1),p(2))-regionLeft); 
                            deltaCRight = abs(regionRight(p(1),p(2))-regionRight); 
                            
                            wLeft = exp(-deltaCLeft/gammaC).*exp(-deltaG/windowDist);
                            wRight = exp(-deltaCRight/gammaC).*exp(-deltaG/windowDist);
                     
                            init_cost = min(abs(regionLeft-regionRight),T);
                            bil = sum(sum(wLeft(:).*wRight(:).*init_cost(:))) / sum(sum(wLeft(:).*wRight(:)));
                            if bil > maxValue
                                maxValue = bil;
                                bestDisplay = s + 1;
                            end
                            
                    end
                end
            end
            disparity(i, j) = abs(bestDisplay);
        end
    end
end

