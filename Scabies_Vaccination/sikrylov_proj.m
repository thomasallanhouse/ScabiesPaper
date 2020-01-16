function [f,errest] = sikrylov_proj(A,b,m,shift)
% SIKRYLOV  Shift-and-invert Arnoldi method for expm(A)*b with projection
% to probability vector.
% Stefan Guettel, 2017.

tol = 0;
if m < 1, % error tolerance
    tol = m;
    m = 50; % maximum number of iter = 30
end

if nargin < 4, shift = 1; end

B = A - shift*speye(size(A));
V = zeros(length(b),m); H = zeros(m);
beta = norm(b);
w = b/beta;
V(:,1) = w;
delay = 2; % for stopping criterion
errest = [];
for c = 1:m,
    w =  B\V(:,c);
    
    % MGS
    for r = 1:c,
        H(r,c) = V(:,r)'*w;
        w = w - V(:,r)*H(r,c);
    end
    
    % CGS
    %H(1:c,c) = V(:,1:c)'*w;
    %w = w - V(:,1:c)*H(1:c,c);

    H(c+1,c) = norm(w);
    V(:,c+1) = w/H(c+1,c);
    
    if tol && c>delay,
        I = eye(c);
        E = expm(H(1:c,1:c)\I); E1 = E(:,1); 
        decay = (beta*exp(shift))*E1;
        errest(c-delay) = norm(decay(end-delay+1:end),inf)/norm(decay,inf);
        if errest(c-delay) < tol, 
            f = V(:,1:c)*decay;
            f(f<0) = 0; f(f>1) = 1; f = f/norm(f,1);
            return
        end
    end

end
warning(['SIKRYLOV_PROJ: Error tolerance not satisfied after ' num2str(c) ' iterations'])
I = eye(c);
E = expm(H(1:c,1:c)\I); E1 = E(:,1); 
decay = (beta*exp(shift))*E1;
f = V(:,1:c)*decay;
f(f<0) = 0; f(f>1) = 1; f = f/norm(f,1);


