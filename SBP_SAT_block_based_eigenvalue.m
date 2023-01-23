%% Initial Parameters
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

Domain = [0,1;0,1];
Bs = 17;
loops = 500;
nplot = 10;

% calculate Eigenvalues with time step evaluation (RK4 = 1) or just the
% FD-operator (RK4 = 0)
RK4 = 0;





%% Treecode setting
% level 2
%treecode = {[0,0],[0,1],[0,2],[0,3],[1,0],[1,1],[1,2],[1,3],[2,0],[2,1],[2,2],[2,3],[3,0],[3,1],[3,2],[3,3]};

% level 2 and 3
%treecode = {[0,0],[0,1,0],[0,1,1],[0,1,2],[0,1,3],[0,2],[0,3],[1,0],[1,1],[1,2],[1,3],[2,0],[2,1],[2,2],[2,3],[3,0],[3,1],[3,2],[3,3]};

% level 3 and 4
%treecode = {[0,0,0],[0,0,1],[0,0,2],[0,0,3,0],[0,0,3,1],[0,0,3,2],[0,0,3,3],[0,1,0],[0,1,1],[0,1,2],[0,1,3],[0,2,0],[0,2,1],[0,2,2],[0,2,3],[0,3,0],[0,3,1],[0,3,2],[0,3,3],[1,0,0],[1,0,1],[1,0,2],[1,0,3],[1,1,0],[1,1,1],[1,1,2],[1,1,3],[1,2,0],[1,2,1],[1,2,2],[1,2,3],[1,3,0],[1,3,1],[1,3,2],[1,3,3],[2,0,0],[2,0,1],[2,0,2],[2,0,3],[2,1,0],[2,1,1],[2,1,2],[2,1,3],[2,2,0],[2,3,0],[3,0,0],[3,1,0],[3,2,0],[3,3,0],[2,2,1],[2,3,1],[3,0,1],[3,1,1],[3,2,1],[3,3,1],[2,2,2],[2,3,2],[3,0,2],[3,1,2],[3,2,2],[3,3,2],[2,2,3],[2,3,3],[3,0,3],[3,1,3],[3,2,3],[3,3,3]};

% level 2, 3 and 4
%treecode = {[0,0,0],[0,0,1],[0,0,2],[0,0,3,0],[0,0,3,1],[0,0,3,2],[0,0,3,3],[0,1,0],[0,1,1],[0,1,2],[0,1,3],[0,2,0],[0,2,1],[0,2,2],[0,2,3],[0,3,0],[0,3,1],[0,3,2],[0,3,3],[1,0],[1,1],[1,2],[1,3],[2,0],[2,1],[2,2],[2,3],[3,0],[3,1],[3,2],[3,3]};

% level 2, 3 and 4 with 2-level-diff over diagonal
%treecode = {[0,0,0],[0,0,1],[0,0,2],[0,0,3,0],[0,0,3,1],[0,0,3,2],[0,0,3,3],[0,1,0],[0,1,1],[0,1,2],[0,1,3],[0,2,0],[0,2,1],[0,2,2],[0,2,3],[0,3],[1,0],[1,1],[1,2],[1,3],[2,0],[2,1],[2,2],[2,3],[3,0],[3,1],[3,2],[3,3]};

% 2,3,4 and 5
treecode = {[0,0,0],[0,0,1],[0,0,2],[0,0,3,0],[0,0,3,1],[0,0,3,2],[0,0,3,3,0],[0,0,3,3,1],[0,0,3,3,2],[0,0,3,3,3],[0,1,0],[0,1,1],[0,1,2,0],[0,1,2,1],[0,1,2,2],[0,1,2,3],[0,1,3],[0,2,0],[0,2,1,0],[0,2,1,1],[0,2,1,2],[0,2,1,3],[0,2,2],[0,2,3],[0,3,0,0],[0,3,0,1],[0,3,0,2],[0,3,0,3],[0,3,1],[0,3,2],[0,3,3],[1,0],[1,1],[1,2],[1,3],[2,0],[2,1],[2,2],[2,3],[3,0],[3,1],[3,2],[3,3]};

% 2,3,4 and 5 y-axis mirrored
%treecode = {[0,0],[0,1],[0,2],[0,3],[1,0],[1,1],[1,2],[1,3],[2,0,0],[2,0,1],[2,0,2],[2,0,3,0],[2,0,3,1],[2,0,3,2],[2,0,3,3],[2,1,0],[2,1,1],[2,1,2,0],[2,1,2,1],[2,1,2,2],[2,1,2,3],[2,1,3],[2,2,0],[2,2,1,0],[2,2,1,1,0],[2,2,1,1,1],[2,2,1,1,2],[2,2,1,1,3],[2,2,1,2],[2,2,1,3],[2,2,2],[2,2,3],[2,3,0,0],[2,3,0,1],[2,3,0,2],[2,3,0,3],[2,3,1],[2,3,2],[2,3,3],[3,0],[3,1],[3,2],[3,3]};

% 2,3,4 and 5 with two refined zones
%treecode = {[0,0,0],[0,0,1],[0,0,2],[0,0,3,0],[0,0,3,1],[0,0,3,2],[0,0,3,3,0],[0,0,3,3,1],[0,0,3,3,2],[0,0,3,3,3],[0,1,0],[0,1,1],[0,1,2,0],[0,1,2,1],[0,1,2,2],[0,1,2,3,0],[0,1,2,3,1],[0,1,2,3,2],[0,1,2,3,3],[0,1,3,0],[0,1,3,1],[0,1,3,2],[0,1,3,3],[0,2,0],[0,2,1,0],[0,2,1,1],[0,2,1,2],[0,2,1,3],[0,2,2],[0,2,3],[0,3,0,0],[0,3,0,1],[0,3,0,2],[0,3,0,3],[0,3,1,0,0],[0,3,1,0,1],[0,3,1,0,2],[0,3,1,0,3],[0,3,1,1],[0,3,1,2],[0,3,1,3],[0,3,2],[0,3,3],[1,0,0],[1,0,1],[1,0,2],[1,0,3],[1,1],[1,2,0],[1,2,1],[1,2,2],[1,2,3],[1,3],[2,0],[2,1],[2,2],[2,3],[3,0],[3,1],[3,2],[3,3]};



nblocks = size(treecode,2);  % number of blocks
maxlevel = max(cellfun(@length, treecode));
minlevel = min(cellfun(@length, treecode));
dx_fine = (Domain(1,2)-Domain(1,1))/2^maxlevel/(Bs-1);  % dx-value of finest block
flow_vel = sqrt(2);
dt = dx_fine/flow_vel;
CFL = flow_vel*dt/dx_fine    % print CFL number

% for further calculations put treecode into an array and fill empty values
% with -1
treecode_array = -ones(nblocks,maxlevel);
for block=1:nblocks
    for l=1:length(treecode{block})
        treecode_array(block,l) = treecode{block}(l);
    end
end


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
Data = generate_Data(Domain,Bs,treecode);

% generate array and block-IDs of neighbours
[neighbours_right_array,neighbours_left_array,neighbours_right_index,neighbours_left_index] = create_neighbours_list(treecode_array,Data);
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
E_0_y = sparse(kron(sparse(eye(Bs)),E_N));
E_N_x = sparse(kron(E_N,sparse(eye(Bs))));
E_N_y = sparse(kron(sparse(eye(Bs)),E_0));



vec_len = Bs^2*nblocks;
Y = sparse(zeros(vec_len));
Y_k = cell(5,1);

% Runge Kutta 4 method params
RK4_a = [0,1/2,1/2,1];


for i=1:vec_len
    e_i = zeros(nblocks*vec_range,1);
    e_i(i) = 1;

    for ki=1:5
        Y_k{ki} = zeros(nblocks*vec_range,1);
    end

    if RK4
        % calculate k-Values for RK4
        for ki = 1:4   
            Y_k{ki+1} = right_hand_side(e_i+dt*RK4_a(ki)*Y_k{ki},Data,treecode_array,neighbours_left_array,neighbours_right_array,neighbours_left_index,neighbours_right_index, ...
            neighbours_top_array, neighbours_bottom_array, neighbours_top_index, neighbours_bottom_index,nblocks,vec_range, Bs, D_x_cell, D_y_cell, H_x_inv_cell, H_y_inv_cell, H_two_blocks, E_0_x, E_0_y, E_N_x, E_N_y, I_F2C, I_C2F,I_single2double);        
        end
    
        % do time step
        Y(:,i) = e_i + 1/6*dt*( Y_k{2} + 2*Y_k{3} + 2*Y_k{4} + Y_k{5});

    else
        Y(:,i) = right_hand_side(e_i,Data,treecode_array,neighbours_left_array,neighbours_right_array,neighbours_left_index,neighbours_right_index, ...
            neighbours_top_array, neighbours_bottom_array, neighbours_top_index, neighbours_bottom_index,nblocks,vec_range, Bs, D_x_cell, D_y_cell, H_x_inv_cell, H_y_inv_cell, H_two_blocks, E_0_x, E_0_y, E_N_x, E_N_y, I_F2C, I_C2F,I_single2double);
    end
    
end

ei = eig(full(Y));

% plot eigenvalues
f = figure;
%ei = eigs(Y,400);
plot(ei,'o')
hold on

if RK4
    theta = 0:0.00001:2*pi;
    plot(cos(theta)+1i*sin(theta))
    title("Eigenvalues of SBP-FD with RK4 Time Step, Level: " +num2str(minlevel)+ " - " +num2str(maxlevel),'interpreter','latex')
else
    title("Eigenvalues of SBP-FD, Level: " +num2str(minlevel)+ " - " +num2str(maxlevel),'interpreter','latex')
end

ylabel("Im",'interpreter','latex')
xlabel("Re",'interpreter','latex')
f.Position = [100 0 1000 1000];
ax = gca; 
ax.FontSize = 20; 


