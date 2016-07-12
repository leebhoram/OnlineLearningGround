function [ V, Zi, Zv, alpha, Vxy, wxz, w_ ] = computeIncMotion( tripleMatchSet, nite )

global f cv cu 
global dt;

V = [];
Z = [];
alpha = 0;

if size(tripleMatchSet,1) < 10
    return;
end

f2 = f*f;

% Initialization using least squares

% Form (A'*A)x = (A'*b)
Idx = [];
B_temp = [];
A_temp = [];
N = size(tripleMatchSet,1);
for m = 1:N 
    u_prev = tripleMatchSet(m,1); v_prev = tripleMatchSet(m,2);
    u_cur = tripleMatchSet(m,3);  v_cur = tripleMatchSet(m,4);
    u_next = tripleMatchSet(m,5); v_next = tripleMatchSet(m,6);
             
    u_prev = u_prev - cu; v_prev = v_prev - cv;
    u_cur = u_cur - cu;    v_cur = v_cur - cv;
    u_next =u_next - cu;   v_next = v_next - cv;
    
    % compensate rotation here
    u_dot = (u_next - u_prev)/(2*dt);
    v_dot = (v_next - v_prev)/(2*dt);
    
    a = (u_cur/f);
    b = u_cur*v_cur/f2;
    c = -(1 + u_cur*u_cur/f2);
    d = v_cur/f;
    
    a_ = d;
    b_ = 1 + d*d;
    c_ = -a*d;
    d_ = -a;
  
    if ( abs(u_cur) > 3) && (abs(v_cur) > 3) 
        A_ = [a b c d; a_ b_ c_ d_];
        B_ = [u_dot/f ; v_dot/f];
        Idx = [Idx; m];
        A_temp = [A_temp; A_];
        B_temp = [B_temp; B_]; 
    end
      
end

r = length(Idx);

A11 = A_temp(:,2:4)'*A_temp(:,2:4);
A11(1,1) = A11(1,1) + 1;
A11(3,3) = A11(3,3) + 1;

A12 = zeros(3,r);
A22 = zeros(r,1);
B = zeros(r+3,1);
B(1:3) = A_temp(:,2:4)'*B_temp;
for n=1:r
    A12(1,n) = A_temp((2*n-1),2)*A_temp((2*n-1),1)+ A_temp((2*n),2).*A_temp((2*n),1);
    A12(2,n) = A_temp((2*n-1),3)*A_temp((2*n-1),1)+ A_temp((2*n),3).*A_temp((2*n),1);
    A12(3,n) = A_temp((2*n-1),4)*A_temp((2*n-1),1)+ A_temp((2*n),4).*A_temp((2*n),1);

    A22(n) = A_temp((2*n-1),1)*A_temp((2*n-1),1) + A_temp((2*n),1).*A_temp((2*n),1);

    B(n+3) = B_temp((2*n-1),1)*A_temp((2*n-1),1) + B_temp((2*n),1).*A_temp((2*n),1);
end

sol_ls = mexSolveArrow3(A11,A12,A22, B);


B_Ax = zeros(2*r,1);
for n=1:r
    B_Ax(2*n-1) = B_temp(2*n-1) -  (A_temp((2*n-1),2:4)*sol_ls(1:3) +  A_temp((2*n-1),1)*sol_ls(3+n));
    B_Ax(2*n  ) = B_temp(2*n  ) -  (A_temp((2*n  ),2:4)*sol_ls(1:3) +  A_temp((2*n  ),1)*sol_ls(3+n));      
end

w_ = sqrt(B_Ax.*B_Ax);    
W = (1./max(w_,1e-5));

% non-linear least squares
sol_nlls = [0; 0; sol_ls];  
J_temp = [zeros(2*r,2) A_temp(:,2:4) A_temp(:,1)];
    
for tt = 1:nite

   
    for n=1:r
        J_temp(2*n-1,1) = -sol_nlls(5+n);   J_temp(2*n-1,6) = -sol_nlls(1) + J_temp(2*n-1,6);
        J_temp(2*n  ,2) = -sol_nlls(5+n);   J_temp(2*n  ,6) = -sol_nlls(2) + J_temp(2*n  ,6);     
    end
 
    AW = [J_temp(:,1).*W J_temp(:,2).*W J_temp(:,3).*W  J_temp(:,4).*W J_temp(:,5).*W]; 
    
    A11 = AW'*J_temp(:,1:5);
    A11(1,1) = A11(1,1);% + 0.1;
    A11(2,2) = A11(2,2);% + 0.1;
    A11(3,3) = A11(3,3);% + 0.1;
    A11(5,5) = A11(5,5);% + 0.1;

    A12 = zeros(5,r);
    A22 = zeros(r,1);
    B = zeros(r+5,1);
    B(1:5) = AW'*B_Ax;
    for n=1:r
        A12(1,n) = J_temp((2*n-1),1)*J_temp((2*n-1),6)*W(2*n-1)+ J_temp((2*n),1).*J_temp((2*n),6)*W(2*n);
        A12(2,n) = J_temp((2*n-1),2)*J_temp((2*n-1),6)*W(2*n-1)+ J_temp((2*n),2).*J_temp((2*n),6)*W(2*n);
        A12(3,n) = J_temp((2*n-1),3)*J_temp((2*n-1),6)*W(2*n-1)+ J_temp((2*n),3).*J_temp((2*n),6)*W(2*n);
        A12(4,n) = J_temp((2*n-1),4)*J_temp((2*n-1),6)*W(2*n-1)+ J_temp((2*n),4).*J_temp((2*n),6)*W(2*n);
        A12(5,n) = J_temp((2*n-1),5)*J_temp((2*n-1),6)*W(2*n-1)+ J_temp((2*n),5).*J_temp((2*n),6)*W(2*n);

        A22(n) = J_temp((2*n-1),6)*J_temp((2*n-1),6)*W(2*n-1) + J_temp((2*n),6).*J_temp((2*n),6)*W(2*n);
        B(n+5) = B_Ax(2*n-1)*J_temp((2*n-1),6)*W(2*n-1) + B_Ax(2*n).*J_temp((2*n),6)*W(2*n);    
    end

    % dsol = solveULArrow(A11,A12,A22, B);
    dsol = mexSolveArrow3(A11,A12,A22, B); 
      
    sol_nlls = sol_nlls + dsol;
   
    
    B_Ax = zeros(2*r,1);
    for n=1:r
        B_Ax(2*n-1) = B_temp(2*n-1) -  (A_temp((2*n-1),2:4)*sol_nlls(3:5) + ( A_temp((2*n-1),1) - sol_nlls(1))*sol_nlls(5+n));
        B_Ax(2*n  ) = B_temp(2*n  ) -  (A_temp((2*n  ),2:4)*sol_nlls(3:5) + ( A_temp((2*n  ),1) - sol_nlls(2))*sol_nlls(5+n));      
    end
    
    if norm(dsol) < 1e-3
        break;
    end
    
    w_ = sqrt(B_Ax.*B_Ax);
    
    W = (1./max(w_,1e-5));
   % W = W/norm(W);
end



alpha = sol_nlls(4);
wxz = sol_nlls([3 5]);

V = 1;
Vxy = sol_nlls(1:2);

Z = [Idx sol_nlls(6:(r+5)) reshape(w_, [2 length(Idx)])' ];
Z = abs(Z(find(sign(Z(:,2)) == sign(median(Z(:,2)))),:));

Zi = Z(:,1);
Zv = Z(:,2);

w_ = Z(:,3:4);

end

