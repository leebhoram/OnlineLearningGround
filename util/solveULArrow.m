function x = solveULArrow_v3(A11,A12,A22,b)
    A11dim = size(A11,1);
    b1 = b(1:A11dim);
    b2 = b((A11dim+1):end);

    u = b2./A22;
    for tt = 1:size(A12,1)
        v(:,tt) = A12(tt,:)'./A22; %#ok<AGROW>
    end

    x = zeros(A11dim+length(u),1);
    x(1:A11dim) = (A11-A12*v)\(b1 - A12*u);
    x((A11dim+1):end) = u-v*x(1:A11dim);
end 


