function U_n = finite_diff_advance(U_n,U_o,K,H,p,tol,fd_F,fd_jac, ...
    bc_L_fun,bc_R_fun,bc_L_jac_fun,bc_R_jac_fun) 

old = U_n;

diff = tol+1;
while diff > tol

    Fout = fd_F(U_n,U_o,K,H,p);

    Jout = fd_jac(U_n,U_o,K,H,p);

    n = size(Fout,1);
    N = size(Fout,2);

    Fval = reshape(Fout,n*N,1);

    % This is slow and should eventually be changed for speed up!
    Jac = sparse(n*(N+2),n*(N+2));
    for j = 1:N % corresponds to node location
        for k = 1:n % corresponds to j-1, j, j+1
            for l = 1:n % corresponds to system equation
                for r = 1:n % corresponds to system variable
                    if Jout{l}{r}{k} == 0
                        Jac(j*n+l,(j-1)*n+(k-1)*n+r) = 0; 
                    else
                        Jac(j*n+l,(j-1)*n+(k-1)*n+r) = Jout{l}{r}{k}(j);
                    end
                end
            end
        end
    end

    bc_L = bc_L_fun(U_n,U_o,K,H,p);
    bc_R = bc_R_fun(U_n,U_o,K,H,p);
    bc_L_jac = bc_L_jac_fun(U_n,U_o,K,H,p);
    bc_R_jac = bc_R_jac_fun(U_n,U_o,K,H,p);

    F = [bc_L;Fval;bc_R];
    Jac(1:size(bc_L_jac,1),1:size(bc_L_jac,2)) = bc_L_jac;
    Jac(end-size(bc_R_jac,1)+1:end,end-size(bc_R_jac,2)+1:end) = bc_R_jac;


    U_n = reshape(reshape(U_n,n*(N+2),1)-Jac\F,n,N+2);
    diff = max(norm(U_n-old));
    old = U_n;

end

