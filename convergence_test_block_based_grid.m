Bs_exp_vec = 4:8;
%% Initial Parameters
Domain = [0,1;0,1];
loops = 4000;
nplot = 1000;
dt = 1/4000;
% characterize the initial data (0 = sin functions, 1 = gauss blob, 2 =
% random function)
init_data = 1;

% Gauss Blob params
sigma = [0.2,0.2];
mu = [0.20,0.25];
%mu = [0.25,0.75];

% sine-function params
load('a_rand_50.mat')
%a = rand(50,1);

% plot parameters
if init_data == 0
    c_range = [-0.2,0.2];
elseif init_data == 1
    c_range = [-0.1,1];
end
    %c_range = [-0.05 0.05];
x_range = [Domain(1,1),Domain(1,2)];
y_range = [Domain(2,1),Domain(2,2)];

% Runge Kutta 4 method params
RK4_a = [0,1/2,1/2,1];


%% Treecode setting
% level 2
%treecode = {[0,0],[0,1],[0,2],[0,3],[1,0],[1,1],[1,2],[1,3],[2,0],[2,1],[2,2],[2,3],[3,0],[3,1],[3,2],[3,3]};

% level 2 and 3
treecode = {[0,0],[0,1,0],[0,1,1],[0,1,2],[0,1,3],[0,2],[0,3],[1,0],[1,1],[1,2],[1,3],[2,0],[2,1],[2,2],[2,3],[3,0],[3,1],[3,2],[3,3]};

% level 3 and 4
%treecode = {[0,0,0],[0,0,1],[0,0,2],[0,0,3,0],[0,0,3,1],[0,0,3,2],[0,0,3,3],[0,1,0],[0,1,1],[0,1,2],[0,1,3],[0,2,0],[0,2,1],[0,2,2],[0,2,3],[0,3,0],[0,3,1],[0,3,2],[0,3,3],[1,0,0],[1,0,1],[1,0,2],[1,0,3],[1,1,0],[1,1,1],[1,1,2],[1,1,3],[1,2,0],[1,2,1],[1,2,2],[1,2,3],[1,3,0],[1,3,1],[1,3,2],[1,3,3],[2,0,0],[2,0,1],[2,0,2],[2,0,3],[2,1,0],[2,1,1],[2,1,2],[2,1,3],[2,2,0],[2,3,0],[3,0,0],[3,1,0],[3,2,0],[3,3,0],[2,2,1],[2,3,1],[3,0,1],[3,1,1],[3,2,1],[3,3,1],[2,2,2],[2,3,2],[3,0,2],[3,1,2],[3,2,2],[3,3,2],[2,2,3],[2,3,3],[3,0,3],[3,1,3],[3,2,3],[3,3,3]};

% level 2, 3 and 4
%treecode = {[0,0,0],[0,0,1],[0,0,2],[0,0,3,0],[0,0,3,1],[0,0,3,2],[0,0,3,3],[0,1,0],[0,1,1],[0,1,2],[0,1,3],[0,2,0],[0,2,1],[0,2,2],[0,2,3],[0,3,0],[0,3,1],[0,3,2],[0,3,3],[1,0],[1,1],[1,2],[1,3],[2,0],[2,1],[2,2],[2,3],[3,0],[3,1],[3,2],[3,3]};

% 2,3,4 and 5
%treecode = {[0,0,0],[0,0,1],[0,0,2],[0,0,3,0],[0,0,3,1],[0,0,3,2],[0,0,3,3,0],[0,0,3,3,1],[0,0,3,3,2],[0,0,3,3,3],[0,1,0],[0,1,1],[0,1,2,0],[0,1,2,1],[0,1,2,2],[0,1,2,3],[0,1,3],[0,2,0],[0,2,1,0],[0,2,1,1],[0,2,1,2],[0,2,1,3],[0,2,2],[0,2,3],[0,3,0,0],[0,3,0,1],[0,3,0,2],[0,3,0,3],[0,3,1],[0,3,2],[0,3,3],[1,0],[1,1],[1,2],[1,3],[2,0],[2,1],[2,2],[2,3],[3,0],[3,1],[3,2],[3,3]};

maxlevel = max(cellfun(@length, treecode));
minlevel = min(cellfun(@length, treecode));
dx_fine = (Domain(1,2)-Domain(1,1))/2^maxlevel/(2^max(Bs_exp_vec));  % dx-value of finest block
flow_vel = sqrt(2);
CFL = flow_vel*dt/dx_fine    % print CFL number

nblocks = size(treecode,2);  % number of blocks
max_Bs = 2^max(Bs_exp_vec)+1;

% for further calculations put treecode into an array and fill empty values
% with -1
treecode_array = -ones(nblocks,maxlevel);
for block=1:nblocks
    for l=1:length(treecode{block})
        treecode_array(block,l) = treecode{block}(l);
    end
end


error_vec_two_norm = zeros(length(Bs_exp_vec),1);
error_vec_infty_norm = zeros(length(Bs_exp_vec),1);
u_Bs = zeros(length(Bs_exp_vec),nblocks*max_Bs^2);

max_array = zeros(size(Bs_exp_vec,2),nblocks)
for exp = 1:size(Bs_exp_vec,2)
    Bs_exp = Bs_exp_vec(exp);
    Bs = 2^Bs_exp+1;
    
    
        %% Operators
    % weight Operators and Derivative Operators (cells for level dependency)
    H_x_inv_cell = cell(maxlevel,1);
    H_y_inv_cell = cell(maxlevel,1);
    % h-weights (without cell step size) for two block interface
    H_two_blocks = zeros(Bs+1);    
    D_x_cell = cell(maxlevel,1);
    D_y_cell = cell(maxlevel,1);
    
    M_cell = cell(maxlevel,1); % For M-norm
    
    % Interpolation Operators
    I_F2C = zeros((Bs+1)/2,Bs);
    I_C2F = zeros(Bs,(Bs+1)/2);
    
    % Projection Operators
    I_single2double = zeros(Bs+1,Bs);   % Projection operator that projects one block onto two blocks
    % (interpolation operator needed)
    E_0_x = zeros(Bs^2,Bs);
    E_0_y = zeros(Bs^2,Bs);
    E_N_x = zeros(Bs*(2*Bs-1),2*Bs-1);
    E_N_y = zeros(Bs*(2*Bs-1),2*Bs-1);
    
    
    %% generate Data and Operators
    %Data = generate_Data_treecode_Wabbit(Domain,Bs,treecode);
    Data = generate_Data(Domain,Bs,treecode);
    
    
    % generate array and block-IDs of neighbours
    [neighbours_right_array,neighbours_left_array,neighbours_right_index,neighbours_left_index] = create_neighbours_list(treecode_array,Data);
    %[neighbours_top_array, neighbours_bottom_array, neighbours_top_index, neighbours_bottom_index] = create_top_neighbours_list(treecode_array,Data);
    [neighbours_top_array, neighbours_bottom_array, neighbours_top_index, neighbours_bottom_index] = create_bottom_neighbours_list(treecode_array,Data);
    
    
    % generate full Data vector for RK-substeps
    vec_range = Bs^2;
    u_k = cell(5,1);
    for ki=1:5
        u_k{ki} = zeros(nblocks*vec_range,1);
    end
     
    % Projection operator that projects one block onto two blocks
    % (interpolation operator needed)
    I_single2double = generate_projection_operator((Bs+1)/2);   
    
    % weight matrix for two block interface
    h_stencil = [17/48,59/48,43/48,49/48];
    H_two_blocks = eye(Bs+1);
    for i=1:4
        H_two_blocks((Bs+1)/2-i+1,(Bs+1)/2-i+1) = h_stencil(i);
        H_two_blocks((Bs+1)/2+i,(Bs+1)/2+i) = h_stencil(i);
    end
    
    % level wise weight matrix and derivataive operators
    for level=1:maxlevel
        [H,D] = generate_SBP_FD_4th(Bs,(Domain(1,2)-Domain(1,1))/(2^(level)*(Bs-1)));
    
        D_x_cell{level} = sparse(kron(D,sparse(eye(Bs))));
        D_y_cell{level} = sparse(kron(sparse(eye(Bs)),D));
        H_x = sparse(kron(H,sparse(eye(Bs))));
        H_y = sparse(kron(sparse(eye(Bs)),H));
        H_x_inv_cell{level} = inv(H_x);
        H_y_inv_cell{level} = inv(H_y);
        % since we only stack finer blocks in the x-axis we have the same size
        % for H and D matrices (x-derivative). But we take 2*Bs-1 for the
        % y-length (see kronecker product)
        M_cell{level} = sparse(kron(H,H));
    
    end
       
    % generate interpolation operators
    [H,D] = generate_SBP_FD_4th((Bs+1)/2,2);
    [H_fine,D_fine] = generate_SBP_FD_4th(Bs,1);
    [I_F2C,I_C2F] = generate_Interpolatin_Operators_4th(Bs,(Bs+1)/2,H_fine,H) ;
    
    % projection operators for SAT-Terms
    E_0 = zeros(Bs,1);
    E_0(1) = 1;
    
    E_N = zeros(Bs,1);
    E_N(end) = 1;
    
    E_0_x = sparse(kron(E_0,sparse(eye(Bs))));
    %E_0_y = sparse(kron(sparse(eye(Bs)),E_0));
    E_0_y = sparse(kron(sparse(eye(Bs)),E_N));
    E_N_x = sparse(kron(E_N,sparse(eye(Bs))));
    %E_N_y = sparse(kron(sparse(eye(Bs)),E_N));
    E_N_y = sparse(kron(sparse(eye(Bs)),E_0));
    
    
    %% create Plot
    clear M
    M(floor(loops/nplot)+1) = struct('cdata',[],'colormap',[]);
    h = figure;
    h.Visible = 'on';
    hold off
    hold on
    colormap turbo;
    colorbar;

    
    for i=1:nblocks
        X = Data.X{i};
        Y = Data.Y{i};
        if init_data == 0   
            Z = sin_cos(X,Y,a);
            max_array(exp,i) = max(max(abs(Z)));
        elseif init_data == 1
            Z = gauss(X,Y,sigma,mu);
        elseif init_data == 2
            Z = rand(Bs,Bs);
        end
        Data.Z{i} = Z;
        imagesc([X(1,1),X(end,end)],[Y(1,1),Y(end,end)],Z);
        xlim(x_range)
        ylim(y_range)
    end
    
    caxis(c_range)
    
    drawnow
    
        
    
    
     
    
    
    
    %% start simulation
    for t=1:loops
        % transform data into vector
        u = zeros(nblocks*vec_range,1);
        for block = 1:nblocks
            u((block-1)*vec_range+1:block*vec_range) = reshape(Data.Z{block},vec_range,1);
        end
        
        % calculate k-Values for RK4
        for ki = 1:4
            %u_k{ki+1} = right_hand_side_treecode_WABBIT(u+dt*RK4_a(ki)*u_k{ki},Data,treecode_array,neighbours_left_array,neighbours_right_array,neighbours_left_index,neighbours_right_index, ...
            %neighbours_top_array, neighbours_bottom_array, neighbours_top_index, neighbours_bottom_index,nblocks,vec_range, Bs, D_x_cell, D_y_cell, H_x_inv_cell, H_y_inv_cell, H_two_blocks, E_0_x, E_0_y, E_N_x, E_N_y, I_F2C, I_C2F,I_single2double);
            u_k{ki+1} = right_hand_side(u+dt*RK4_a(ki)*u_k{ki},Data,treecode_array,neighbours_left_array,neighbours_right_array,neighbours_left_index,neighbours_right_index, ...
        neighbours_top_array, neighbours_bottom_array, neighbours_top_index, neighbours_bottom_index,nblocks,vec_range, Bs, D_x_cell, D_y_cell, H_x_inv_cell, H_y_inv_cell, H_two_blocks, E_0_x, E_0_y, E_N_x, E_N_y, I_F2C, I_C2F,I_single2double);
    
        end
        
        % do time step
        u = u+1/6*dt*(u_k{2} + 2*u_k{3} + 2*u_k{4} + u_k{5});
        for block = 1:nblocks
            Data.Z{block} = reshape(u((block-1)*vec_range+1:block*vec_range),Bs,Bs);
        end
    
    
        %% plot
        
        if mod(t,nplot) == 0
            hold off
            hold on 
            h.Visible = 'on';
            colormap turbo;
            colorbar;
            
            for block=1:nblocks
                imagesc([Data.X{block}(1,1),Data.X{block}(end,end)],[Data.Y{block}(1,1),Data.Y{block}(end,end)],Data.Z{block})
                xlim(x_range)
                ylim(y_range)
            end
    
            caxis(c_range)       
            drawnow
            M(floor(t/nplot)+1) = getframe();
        end
        
    
    
    
    end
    
    
    
    
    solution = zeros(nblocks*vec_range,1);
    for block = 1:nblocks
        if init_data == 0
            solution((block-1)*vec_range+1:block*vec_range) = reshape(sin_cos(Data.X{block}-dt*loops,Data.Y{block},a),vec_range,1);
        else
            solution((block-1)*vec_range+1:block*vec_range) = reshape(gauss(Data.X{block}-dt*loops,Data.Y{block},sigma,mu),vec_range,1);
        end
    end
    
    u_Bs(exp,1:nblocks*vec_range) = u;
    error_vec_two_norm(exp) = norm(u-solution)/norm(solution);
    error_vec_infty_norm(exp) = max(abs(u-solution));


end
 

f = figure;
loglog(Bs_exp_vec,error_vec_two_norm)
if init_data == 0
    title("2-Norm of Error, Level: " +num2str(minlevel)+ " - " +num2str(maxlevel) +", 50 Modes, Steps = "+num2str(loops),'interpreter','latex')
elseif init_data == 1
    title("2-Norm of Error, Level: " +num2str(minlevel)+ " - " +num2str(maxlevel) +", Gauss Blob, Steps = "+num2str(loops),'interpreter','latex')
end
xticks(Bs_exp_vec)
xlabel("Bs exponent to base 2",'interpreter','latex')
ylabel("$\frac{||u-\bar{u}||_2}{||\bar{u}||_2}$",'interpreter','latex')

f.Position = [100 0 1000 1000];
ax = gca; 
ax.FontSize = 20; 


save_cell{1} = u_Bs;
save_cell{2} = error_vec_two_norm;
save_cell{3} = error_vec_infty_norm;

if init_data == 1
    save("error_infty_exp_"+num2str(Bs_exp_vec(1))+"_"+num2str(Bs_exp_vec(end))+"_level_"+num2str(minlevel)+"_"+num2str(maxlevel)+"_gauss_loops="+num2str(loops)+".mat",'save_cell');
elseif init_data == 0
    save("error_infty_exp_"+num2str(Bs_exp_vec(1))+"_"+num2str(Bs_exp_vec(end))+"_level_"+num2str(minlevel)+"_"+num2str(maxlevel)+"_sin_loops="+num2str(loops)+".mat",'save_cell');
end





