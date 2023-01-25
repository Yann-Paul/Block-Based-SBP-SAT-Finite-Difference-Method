Bs_exp_vec = 5:8;

Domain = [0,1;0,1];
error_vec = zeros(size(Bs_exp_vec,2),1);

c_range = [-0.3,0.3];
x_range = [Domain(1,1),Domain(1,2)];
y_range = [Domain(2,1),Domain(2,2)];

% characterize the initial data (0 = sin functions, 1 = gauss blob, 2 =
% random function)
init_data = 1;

sigma = [0.05 0.05];
mu = [0.5,0.5];
dt = 1/2000;
loops = 2000;
nplot = 200;
dx = (Domain(1,2)-Domain(1,1))/(2*(2*2^Bs_exp_vec(end)));
CFL = dt/dx

load('a_rand_50.mat')
fourier = 1;

for exp = 1:size(Bs_exp_vec,2)
    Bs_exp = Bs_exp_vec(exp);
    Bs = 2^Bs_exp+1;

    
    [H_l,D_l] = generate_SBP_FD_4th(Bs,(Domain(1,2)-Domain(1,1))/(2*(Bs-1)));
    [H_r,D_r] = generate_SBP_FD_4th(2*Bs-1,(Domain(1,2)-Domain(1,1))/(2*(2*Bs-2)));
    
    [I_F2C,I_C2F] = generate_Interpolatin_Operators_4th(2*Bs-1,Bs,H_r,H_l) ;
    
    
    %% Advection equation
    % set equidistant grid and gauss blob
    x_l = linspace(Domain(1,1),Domain(1,1)+1/2*(Domain(1,2)-Domain(1,1)),Bs);
    x_r = linspace(Domain(1,1)+1/2*(Domain(1,2)-Domain(1,1)),Domain(1,2),2*Bs-1);
    y_l = linspace(Domain(2,1),Domain(2,2),Bs);
    y_r = linspace(Domain(2,1),Domain(2,2),2*Bs-1);
    [X_l,Y_l] = meshgrid(x_l,y_l);
    [X_r,Y_r] = meshgrid(x_r,y_r);

    if fourier == 1
        Z_l = sin_cos(X_l,Y_l,a);
        Z_r = sin_cos(X_r,Y_r,a);
    elseif fourier == 0
        Z_l = gauss(X_l,Y_l,sigma,mu);
        Z_r = gauss(X_r,Y_r,sigma,mu);
    end
     

    
    [numRows_l,numCols_l] = size(Z_l);
    vec_len_l = numRows_l*numCols_l;
    u = reshape(Z_l,vec_len_l,1);
    
    [numRows_r,numCols_r] = size(Z_r);
    vec_len_r = numRows_r*numCols_r;
    v = reshape(Z_r,vec_len_r,1);
    
    D_L = sparse(kron(D_l,sparse(eye(numRows_l))));
    H_L = sparse(kron(H_l,sparse(eye(numRows_l))));
    H_L_inv = inv(H_L);
    e_0_l = zeros(numCols_l,1);
    e_0_l(1) = 1;
    e_0_L = sparse(kron(e_0_l,sparse(eye(numRows_l))));
    e_n_l = zeros(numCols_l,1);
    e_n_l(end) = 1;
    e_n_L = sparse(kron(e_n_l,sparse(eye(numRows_l))));
    
    D_R = sparse(kron(D_r,sparse(eye(numRows_r))));
    H_R = sparse(kron(H_r,sparse(eye(numRows_r))));
    H_R_inv = inv(H_R);
    e_0_r = zeros(numCols_r,1);
    e_0_r(1) = 1;
    e_0_R = sparse(kron(e_0_r,sparse(eye(numRows_r))));
    e_n_r = zeros(numCols_r,1);
    e_n_r(end) = 1;
    e_n_R = sparse(kron(e_n_r,sparse(eye(numRows_r))));
    
    M_L = sparse(kron(H_l,H_l));
    M_R = sparse(kron(H_r,H_r));

    
    % plot the original gauss blob first
    clear M
    M(floor(loops/nplot)+1) = struct('cdata',[],'colormap',[]);
    h = figure;
    h.Visible = 'on';
    colormap turbo;
    colorbar;
    pcolor(X_l,Y_l,Z_l)
    caxis(c_range);

    %M(1) = getframe();

    u_norm = zeros(loops+1,1);
    u_norm_0 = u'*M_L*u;

    for t=1:loops
    
        u_0 = u(1:numRows_l);
        u_N = u(vec_len_l-numRows_l+1:vec_len_l);
    
        v_0 = v(1:numRows_r);
        v_N = v(vec_len_r-numRows_r+1:vec_len_r);
        
        k1_L = (-D_L*u+1/2.*H_L_inv*(e_n_L*(u_N-I_F2C*v_0))-1/2.*H_L_inv*(e_0_L*(u_0-I_F2C*v_N)));
        k1_R = (-D_R*v+1/2.*H_R_inv*(e_n_R*(v_N-I_C2F*u_0))-1/2.*H_R_inv*(e_0_R*(v_0-I_C2F*u_N)));
        %k1_L = -D_L*u+1/2.*H_L_inv*(e_n_L*(u_N-I_F2C*v_0));
        %k1_R = -D_R*v-1/2.*H_R_inv*(e_0_R*(v_0-I_C2F*u_N));
    
        u_new = u + 1/2 *dt* k1_L;
        u_0_new = u_new(1:numRows_l);
        u_N_new = u_new(vec_len_l-numRows_l+1:vec_len_l);
    
        v_new = v + 1/2 *dt* k1_R;
        v_0_new = v_new(1:numRows_r);
        v_N_new = v_new(vec_len_r-numRows_r+1:vec_len_r);
    
        k2_L = (-D_L*(u_new)+1/2.*H_L_inv*(e_n_L*(u_N_new-I_F2C*v_0_new))-1/2.*H_L_inv*(e_0_L*(u_0_new-I_F2C*v_N_new)));
        k2_R = (-D_R*(v_new)+1/2.*H_R_inv*(e_n_R*(v_N_new-I_C2F*u_0_new))-1/2.*H_R_inv*(e_0_R*(v_0_new-I_C2F*u_N_new)));
        %k2_L = -D_L*u_new+1/2.*H_L_inv*(e_n_L*(u_N_new-I_F2C*v_0_new));
        %k2_R = -D_R*v_new-1/2.*H_R_inv*(e_0_R*(v_0_new-I_C2F*u_N_new));
    
        u_new = u + 1/2 *dt* k2_L;
        u_0_new = u_new(1:numRows_l);
        u_N_new = u_new(vec_len_l-numRows_l+1:vec_len_l);
    
        v_new = v + 1/2 *dt* k2_R;
        v_0_new = v_new(1:numRows_r);
        v_N_new = v_new(vec_len_r-numRows_r+1:vec_len_r);
    
        k3_L = (-D_L*(u_new)+1/2.*H_L_inv*(e_n_L*(u_N_new-I_F2C*v_0_new))-1/2.*H_L_inv*(e_0_L*(u_0_new-I_F2C*v_N_new)));
        k3_R = (-D_R*(v_new)+1/2.*H_R_inv*(e_n_R*(v_N_new-I_C2F*u_0_new))-1/2.*H_R_inv*(e_0_R*(v_0_new-I_C2F*u_N_new)));
        %k3_L = -D_L*u_new+1/2.*H_L_inv*(e_n_L*(u_N_new-I_F2C*v_0_new));
        %k3_R = -D_R*v_new-1/2.*H_R_inv*(e_0_R*(v_0_new-I_C2F*u_N_new));
    
        u_new = u + dt* k3_L;
        u_0_new = u_new(1:numRows_l);
        u_N_new = u_new(vec_len_l-numRows_l+1:vec_len_l);
    
        v_new = v + dt*k3_R;
        v_0_new = v_new(1:numRows_r);
        v_N_new = v_new(vec_len_r-numRows_r+1:vec_len_r);
        
        k4_L = (-D_L*(u_new)+1/2.*H_L_inv*(e_n_L*(u_N_new-I_F2C*v_0_new))-1/2.*H_L_inv*(e_0_L*(u_0_new-I_F2C*v_N_new)));
        k4_R = (-D_R*(v_new)+1/2.*H_R_inv*(e_n_R*(v_N_new-I_C2F*u_0_new))-1/2.*H_R_inv*(e_0_R*(v_0_new-I_C2F*u_N_new)));
        %k4_L = -D_L*u_new+1/2.*H_L_inv*(e_n_L*(u_N_new-I_F2C*v_0_new));
        %k4_R = -D_R*v_new-1/2.*H_R_inv*(e_0_R*(v_0_new-I_C2F*u_N_new));
    
        u = u + 1/6*dt*(k1_L+ 2*k2_L+ 2*k3_L+ k4_L);
        v = v + 1/6*dt*(k1_R+ 2*k2_R+ 2*k3_R+ k4_R);
    
    %     k1_R = dt*(-D_R*v+1/2.*H_R_inv*(e_n_R*(v_N-I_C2F*u_0))-1/2.*H_R_inv*(e_0_R*(v_0-I_C2F*u_N)));
    %     k2_R = dt*(-D_R*(v+1/2*k1_R)+1/2.*H_R_inv*(e_n_R*(v_N-I_C2F*u_0))-1/2.*H_R_inv*(e_0_R*(v_0-I_C2F*u_N)));
    %     k3_R = dt*(-D_R*(v+1/2*k2_R)+1/2.*H_R_inv*(e_n_R*(v_N-I_C2F*u_0))-1/2.*H_R_inv*(e_0_R*(v_0-I_C2F*u_N)));
    %     k4_R = dt*(-D_R*(v+k3_R)+1/2.*H_R_inv*(e_n_R*(v_N-I_C2F*u_0))-1/2.*H_R_inv*(e_0_R*(v_0-I_C2F*u_N)));
    % 
    %     v = v + 1/6*(k1_R+ 2*k2_R+ 2*k3_R+ k4_R);
    
    
        Z_l = reshape(u,numRows_l,numCols_l);
    
        Z_r = reshape(v,numRows_r,numCols_r);
        
        if mod(t,nplot) == 0
    
            hold off
            h.Visible = 'on';
            colormap turbo;
            colorbar;
            p = pcolor(X_l,Y_l,Z_l);
            set(p, 'EdgeAlpha', '0.1');
            hold on
            p = pcolor(X_r,Y_r,Z_r);
            set(p, 'EdgeAlpha', '0.1');
            xlim(x_range)
            ylim(y_range)
            caxis(c_range);
            colorbar
            %view(0,90);
            drawnow
            M(floor(t/nplot)+1) = getframe();
    
        end

%         u_norm(t+1) = u'*M_L*u-u_norm_0;
    
    end



    mu_sol = [0.5+dt*loops,0.5];

    if fourier == 1
        solution_l = sin_cos(X_l-dt*loops,Y_l,a);
        solution_r = sin_cos(X_r-dt*loops,Y_r,a);
    else
        solution_l = gauss(X_l,Y_l,sigma,mu_sol);
        solution_r = gauss(X_r,Y_r,sigma,mu_sol);
    end

    
    u_v = [u;v];




    solution_l = reshape(solution_l,[],1);
    solution_r = reshape(solution_r,[],1);
    solution = [solution_l;solution_r];
    error_vec(exp) = norm(u_v-solution)/norm(solution);

    
    
    
    
    
    %movie(M);
    
%     
%     figure
%     plot(0:1:loops,u_norm)
%     title("M-norm of solution")
%     xlabel("dt = "+num2str(dt))
    
end

