function [x, num_its, nrm_res] = crs_cg( row_ptr, col_ind, val, y, x, maxit, tol)
% Conjugate Gradient Method for symmetric CRS sparse matrices in NumGeom.
% Algorithm from Trefethen and Bau p. 294.

% %#codegen -args {coder.typeof(int32(0), [Inf,1]), coder.typeof(int32(0), [Inf, 1]), coder.typeof(double(0), [Inf, 1]),
% %#codegen coder.typeof(double(0), [Inf, 1]), coder.typeof(double(0), [Inf, 1]), int32(0), double(0)}

for i=int32(1):length(x)
    x(i) = 0;
end

res = y;

nrm_res = sqrt(res'*res);

search_dir = res;

num_its = int32(0);

Amult_searchdir = nullcopy(res);

while num_its<maxit
    
    Amult_searchdir = crs_prodSymmatVec(row_ptr, col_ind, val, search_dir, Amult_searchdir);
    
    alpha = res'*res/(search_dir'*Amult_searchdir);
    
    x = x + alpha*search_dir;
    
    newres = res - alpha*Amult_searchdir;
    
    beta = newres'*newres/(res'*res);
    
    search_dir = newres + beta*search_dir;
    
    res = newres;
    
    nrm_res = sqrt(res'*res);
    
    num_its = num_its + 1;
    
    if nrm_res<tol
        return;
    end
end
