%% Initial Parameters
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

Domain = [0,1;0,1];
Bs = 17;
loops = 300000;
nplot = 1e2;
dt = 1/4000;
% characterize the initial data (0 = sin functions, 1 = gauss blob, 2 =
% random function)
init_data = 1;

% Gauss Blob params
sigma = [0.2,0.2];
mu = [0.25,0.25]; 

% sine-function params
load('a_rand_50.mat')
%a = rand(30,1);

% plot parameters
if init_data == 1
    c_range = [-0.1,1.1];
elseif init_data == 0
    c_range = [-0.3,0.5];
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


%% create Plot
clear M
M(floor(loops/nplot)+1) = struct('cdata',[],'colormap',[]);
h = figure;
h.Visible = 'on';
hold on
colormap turbo;
colorbar;

for i=1:nblocks
    X = Data.X{i};
    Y = Data.Y{i};
    if init_data == 0   
        Z = sin_cos(X,Y,a);
    elseif init_data == 1
        Z = gauss(X,Y,sigma,mu);
    elseif init_data == 2
        Z = rand(Bs,Bs);
    end
    Data.Z{i} = Z;
    imagesc([X(1,1),X(end,end)],[Y(1,1),Y(end,end)],Z,AlphaData = 0.95);
    %p = pcolor(X,Y,Z);
    %set(p, 'EdgeAlpha', '0.5');
    xlim(x_range)
    ylim(y_range)
end

if init_data == 1
    title("Gauss Blob Initialization, Level: " +num2str(minlevel)+ " - " +num2str(maxlevel) + ", step 0",'interpreter','latex')
elseif init_data == 0
    title("50 Modes Initialization, Level: " +num2str(minlevel)+ " - " +num2str(maxlevel) + ", step 0",'interpreter','latex')
end
ylabel("y",'interpreter','latex')
xlabel("x",'interpreter','latex')
h.Position = [100 0 1000 1000];
ax = gca; 
ax.FontSize = 20; 

caxis(c_range)

M(1) = getframe();




Z_norm = zeros(loops+1,1);
Z_norm_0 = 0;
for block=1:nblocks
    u_block = reshape(Data.Z{block},vec_range,1);
    level = Data.level{block};
    Z_norm_0 = Z_norm_0 + u_block'*M_cell{level}*u_block;
end


%% start simulation
for t=1:loops
    % transform data into vector
    u = zeros(nblocks*vec_range,1);
    for block = 1:nblocks
        u((block-1)*vec_range+1:block*vec_range) = reshape(Data.Z{block},vec_range,1);
    end
    
    % calculate k-Values for RK4
    for ki = 1:4
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
        
        X = zeros(nblocks,2);
        Y = zeros(nblocks,2);
        Z = zeros(nblocks,Bs,Bs);
        for block=1:nblocks
            %X(block,:) = [Data.X{block}(1,1),Data.X{block}(end,end)];
            %Y(block,:) = [Data.Y{block}(1,1),Data.Y{block}(end,end)];
            %Z(block,:,:) = Data.Z{block};
            imagesc([Data.X{block}(1,1),Data.X{block}(end,end)],[Data.Y{block}(1,1),Data.Y{block}(end,end)],Data.Z{block},AlphaData = 0.8)
            xlim(x_range)
            ylim(y_range)
        end

        %imagesc(X,Y,Z,AlphaData = 0.8)
%         for block=4:7
%             imagesc([Data.X{block}(1,1),Data.X{block}(end,end)],[Data.Y{block}(1,1),Data.Y{block}(end,end)],Data.Z{block},AlphaData = 0.8)
%             xlim([0.125,0.25])
%             ylim([0.125,0.25])
%         end


        caxis(c_range)       
        drawnow
        M(floor(t/nplot)+1) = getframe();
    end

    for block=1:nblocks
        u_block = reshape(Data.Z{block},vec_range,1);
        level = Data.level{block};
        u_norm = u_block'*M_cell{level}*u_block;
        Z_norm(t+1) = Z_norm(t+1) + u_block'*M_cell{level}*u_block;
    end
    Z_norm(t+1) = Z_norm(t+1)-Z_norm_0;


end


%movie(M,1, 10)


f = figure
plot(0:1:loops,Z_norm,'k.')
title("M-norm, SBP-FD, dt = " +num2str(dt,'%.1e') +", Level: 2-5",'interpreter','latex')
xlabel("Steps",'interpreter','latex')
ylabel("M-norm",'interpreter','latex')
f.Position = [100 0 1000 1000];
ax = gca; 
ax.FontSize = 20; 

