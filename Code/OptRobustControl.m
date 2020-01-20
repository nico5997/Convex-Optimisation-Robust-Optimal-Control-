function [F, phi_Q, g, C, d_c, x, gam] = OptRobustControl(A,B,C,D,d,Q,R,S,T,L,N,sigma,u0,x0)

n = size(A);
n = n(1);
m = size(B);
m = m(2);
p = length(d);

A_L = zeros(n+1,L*(n + 1)); 
B_L = [zeros(n,L*m); ones(1,L*m)];

% (1) obtention of the differents indexes i of alpha1, alpha12 , etc by taking a 
% random normally distributed term with relative standard deviation
% (2) transormation of the matrix A -> Phi to discretise the system
for i = 1:L
    for j = 1:n
        for p = 1:n
                A_L(j,3*(i-1)+p) = normrnd(A(j,p),sigma*abs(A(j,p))); 
        end
    end
    A_L(1:n,1 + (i-1)*(n+1):2 + (i-1)*(n+1)) = expm(A_L(1:n,1 + (i-1)*(n+1):2 + (i-1)*(n+1))*T); 
end

% (1) obtention of the differents indexes i of beta by taking a 
% random normally distributed term with relative standard deviation
% (2) transormation of the matrix B -> Gamma to discretise the system
for i = 1:L
     for j = 1:n
        for p = 1:m
                B_L(j,p-m+i*m) = normrnd(B(j,p),sigma*abs(B(j,p))); 
        end
     end
    B_L(1:n,1 + m*(i-1):i*m) = inv(A_L(1:n,1 + (i-1)*(n+1):2 + (i-1)*(n+1)))*(A_L(1:n,1 + (i-1)*(n+1):2 + (i-1)*(n+1))*T - eye(n,n))*B_L(1:n,i*m); 
end

n = n + 1;

% concatenation by Block diagonal of the Matrix A
A_Blk = [];
for i = 1:L
    A_Blk = blkdiag(A_Blk, A_L(:,1 + n*(i-1):i*n)); 
end

% concatenation by Block diagonal of the Matrix C
C_Blk = [];
for i = 1:L
    C_Blk = blkdiag(C_Blk, C); 
end

% concatenation of Matrix B
B_Blk = [];
for i = 1:L
    B_Blk = [B_Blk; B_L(:,i)]; 
end

% concatenation of Matrix D
D_tmp = [];
D_Blk = [];
for i = 1:L
    D_tmp = D;
    D_Blk = [D_Blk; D_tmp]; 
end

% concatenation of vector d
d_tmp = []; 
d_Blk = [];
for i = 1:L
    d_tmp = d;
    d_Blk = [d_Blk; d_tmp]; 
end

% Quadraric inegality constraint Matrix Qi concatenation
E_tmp = [];
Q_Blk = [];
phi_Q = [];

for i = 1:L
    for j = 1:N+1
        E_tmp = [zeros(n+m,n*(i-1)), [eye(n,n); zeros(1,n)], zeros(n+m,n*L-n*(i-1)-n),[zeros(n,m); eye(m,m)]];
        Q_bar = E_tmp'*[Q S; S' R]*E_tmp;
        Q_Blk = blkdiag(Q_Blk, Q_bar);
    end
    phi_Q = [phi_Q Q_Blk];
    Q_Blk = [];
end

% Linaer equality constraint Matrix F concatenation
F = [];
for i = 1:N
    F = blkdiag(F, [-A_Blk, -B_Blk]); 
end

tmp1 = size(F);
tmp2 = size(A_Blk);

F = [zeros(tmp2(1),tmp1(2)+tmp2(2)+1); F, zeros(tmp1(1),tmp2(2)+1)];

tmp1 = size(F);
tmp3 = 0;
for i = 1:(tmp1(2)-N*m-1)/tmp2(2)
    for j = 1:tmp1(1)/tmp2(2)
        if (j == i )
        F((j-1)*tmp2(1)+1:j*tmp2(1),(i-1)*(tmp2(2)+m)+1:i*tmp2(2)+(i-1)*m) = eye(tmp2);
        tmp3 = tmp3 + 1;
        end
    end
end

% Linaer inequality constraint Matrix C concatenation
C = [];
for i = 1:N+1
    C = blkdiag(C,[C_Blk, D_Blk]);
end

% Linaer equality constraint vector g concatenation
g = [];
for i = 1:L
g = [g; x0];
end
g = [g; zeros(N*tmp2(1),1)];
g(n*L+1) = 4-1.82;

% Linaer inequality constraint vector d concatenation
d_tmp = [];
d_c = [];
for i = 1:N+1
    d_tmp = d_Blk;
    d_c = [d_c; d_tmp];
end

% CVX optimisation
cvx_begin
    variables gam(1) x((N+1)*(L*n+m));
    minimise(gam)
    subject to
        F*x == g;
        C*x <= d_c;
        x(n*L+1) == u0; % initial input condition added to constraints based on paper graph
        for i = 1:L
           x'*phi_Q(:,(i-1)*length(x)+1:i*length(x))*x <= gam;
        end         
cvx_end

end