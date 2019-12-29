function Error = gs_errfunction(P0, Xobs)
    
    H = reshape(P0(1:9), 3, 3);
    x1 = reshape(Xobs(1:(length(Xobs)/2)), 2, length(Xobs)/4);
    x2 = reshape(Xobs(((length(Xobs)/2)+1):end), 2, length(Xobs)/4);

    Error= [x2 - euclid(H*[x1; ones(1, length(Xobs)/4)]), x1 - euclid(H\[x2; ones(1, length(Xobs)/4)])];
   
end