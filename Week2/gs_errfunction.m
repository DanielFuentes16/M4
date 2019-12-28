function E = gs_errfunction(P0, Xobs)
    
    H  = reshape(P0(1:9),3,3);
    
    Xobs = reshape (Xobs,2,[]);
    x = [Xobs(:,1:size(Xbos,2)/2); ones(1,size(Xbos,2)/2)];
    xp = [Xobs(:, size(Xbos,2)/2+1:end); ones(1,size(Xbos,2)/2)];
    
    xh = P0(9+1:end);
    xhe = [xh(1: length(xh)/2),xh(length(xh)/2+1:end),ones(length(xh)/2,1)];
    xhp = H*xhe';
    
    E = sum((x-xhe).^2)+sum((xp-xhp').^2);
   
end