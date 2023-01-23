function u_new = right_hand_side(u,Data,treecode_array,neighbours_left_array,neighbours_right_array,neighbours_left_index,neighbours_right_index, ...
        neighbours_top_array, neighbours_bottom_array, neighbours_top_index, neighbours_bottom_index,nblocks,vec_range, Bs, D_x_cell, D_y_cell, H_x_inv_cell, H_y_inv_cell, H_two_blocks, E_0_x, E_0_y, E_N_x, E_N_y, I_F2C, I_C2F,I_single2double)
    
    u_new = zeros(size(u));
    for block = 1:nblocks
        local_treecode = treecode_array(block,:);
        level = Data.level{block};
    
        Z = reshape(u((block-1)*vec_range+1:block*vec_range),Bs,Bs);
        u_block = u((block-1)*vec_range+1:block*vec_range);
        u_0 = Z(:,1);
        u_N = Z(:,end);

        u__0 = Z(end ,:)';
        u__N = Z(1,:)';
    
        left_neighbours = neighbours_left_array(block,:,:);
        right_neighbours = neighbours_right_array(block,:,:);
        left_index = neighbours_left_index(block,:);
        right_index = neighbours_right_index(block,:);

        bottom_neighbours = neighbours_bottom_array(block,:,:);
        top_neighbours = neighbours_top_array(block,:,:);
        bottom_index = neighbours_bottom_index(block,:);
        top_index = neighbours_top_index(block,:);
    
        %% generate SAT_L Term
        if Data.level{left_index(1)} == level % if left neighbour has same level
            Z_left = reshape(u((left_index(1)-1)*vec_range+1:left_index(1)*vec_range),Bs,Bs);
            v_N = Z_left(:,end);
            SAT_L = -1/2.*H_x_inv_cell{level}*(E_0_x*(u_0-v_N));
            
        elseif Data.level{left_index(1)} > level % if left neighbour has higher level (two neighbours)
            Z_left_1 = reshape(u((left_index(1)-1)*vec_range+1:left_index(1)*vec_range),Bs,Bs);
            v_N_1 = Z_left_1(:,end);
            Z_left_2 = reshape(u((left_index(2)-1)*vec_range+1:left_index(2)*vec_range),Bs,Bs);
            v_N_2 = Z_left_2(:,end);
            if sum(left_neighbours(:,:,1))<sum(left_neighbours(:,:,2))
                v_N = [I_F2C*v_N_1;I_F2C*v_N_2];
            else
                v_N = [I_F2C*v_N_2;I_F2C*v_N_1];
            end
            SAT_L = -1/2.*H_x_inv_cell{level}*(E_0_x*(u_0-I_single2double'*H_two_blocks*v_N));
            
        else % if left neighbour has lower level
            Z_left = reshape(u((left_index(1)-1)*vec_range+1:left_index(1)*vec_range),Bs,Bs);
            v_N = Z_left(:,end);
            v_n_two_blocks = I_single2double*v_N;
            if local_treecode(level) <= 1 % if this block is upper part
                SAT_L = -1/2.*H_x_inv_cell{level}*(E_0_x*(u_0-I_C2F*v_n_two_blocks(1:(Bs+1)/2)));
            else % if this block is lower part
                SAT_L = -1/2.*H_x_inv_cell{level}*(E_0_x*(u_0-I_C2F*v_n_two_blocks((Bs+1)/2+1:end)));
            end
            
        end
        
        %% generate SAT_R Term
        if Data.level{right_index(1)} == level % if right neighbour has same level
            Z_right = reshape(u((right_index(1)-1)*vec_range+1:right_index(1)*vec_range),Bs,Bs);
            w_0 = Z_right(:,1);
            SAT_R = 1/2.*H_x_inv_cell{level}*(E_N_x*(u_N-w_0));
        elseif Data.level{right_index(1)} > level % if right neighbour has higher level (two neighbours)
            Z_right_1 = reshape(u((right_index(1)-1)*vec_range+1:right_index(1)*vec_range),Bs,Bs);
            w_0_1 = Z_right_1(:,1);
            Z_right_2 = reshape(u((right_index(2)-1)*vec_range+1:right_index(2)*vec_range),Bs,Bs);
            w_0_2 = Z_right_2(:,1);
            if sum(right_neighbours(:,:,1))<sum(right_neighbours(:,:,2))
                w_0 = [I_F2C*w_0_1;I_F2C*w_0_2];
            else
                w_0 = [I_F2C*w_0_2;I_F2C*w_0_1];
            end
            SAT_R = 1/2.*H_x_inv_cell{level}*(E_N_x*(u_N-I_single2double'*H_two_blocks*w_0));
        else % if right neighbour has lower level
            Z_right = reshape(u((right_index(1)-1)*vec_range+1:right_index(1)*vec_range),Bs,Bs);
            w_0 = Z_right(:,1);
            w_0_two_blocks = I_single2double*w_0;
            if local_treecode(level) <= 1 % if this block is upper part
                SAT_R = 1/2.*H_x_inv_cell{level}*(E_N_x*(u_N-I_C2F*w_0_two_blocks(1:(Bs+1)/2)));
            else % if this block is lower part
                SAT_R = 1/2.*H_x_inv_cell{level}*(E_N_x*(u_N-I_C2F*w_0_two_blocks((Bs+1)/2+1:end)));
            end
            
        end

        %% generate SAT_T Term  Top

        if Data.level{top_index(1)} == level % if left neighbour has same level
            Z_top = reshape(u((top_index(1)-1)*vec_range+1:top_index(1)*vec_range),Bs,Bs);
            v__N = Z_top(1,:)';
            SAT_T = -1/2.*H_y_inv_cell{level}*(E_0_y*(u__0-v__N));
            
        elseif Data.level{top_index(1)} > level % if top neighbour has higher level (two neighbours)
            Z_top_1 = reshape(u((top_index(1)-1)*vec_range+1:top_index(1)*vec_range),Bs,Bs);
            v__N_1 = Z_top_1(1,:)';
            Z_top_2 = reshape(u((top_index(2)-1)*vec_range+1:top_index(2)*vec_range),Bs,Bs);
            v__N_2 = Z_top_2(1,:)';
            if sum(top_neighbours(:,:,1))<sum(top_neighbours(:,:,2))
                v__N = [I_F2C*v__N_1;I_F2C*v__N_2];
            else
                v__N = [I_F2C*v__N_2;I_F2C*v__N_1];
            end
            SAT_T = -1/2.*H_y_inv_cell{level}*(E_0_y*(u__0-I_single2double'*H_two_blocks*v__N));
            
        else % if top neighbour has lower level
            Z_top = reshape(u((top_index(1)-1)*vec_range+1:top_index(1)*vec_range),Bs,Bs);
            v__N = Z_top(1,:)';
            v__N_two_blocks = I_single2double*v__N;
            if mod(local_treecode(level),2) == 0 % if this block is left part
                SAT_T = -1/2.*H_y_inv_cell{level}*(E_0_y*(u__0-I_C2F*v__N_two_blocks(1:(Bs+1)/2)));
            else % if this block is right part
                SAT_T = -1/2.*H_y_inv_cell{level}*(E_0_y*(u__0-I_C2F*v__N_two_blocks((Bs+1)/2+1:end)));
            end
            
        end

        %% generate SAT_B Term  Bottom
        if Data.level{bottom_index(1)} == level % if bottom neighbour has same level
            Z_bottom = reshape(u((bottom_index(1)-1)*vec_range+1:bottom_index(1)*vec_range),Bs,Bs);
            w__0 = Z_bottom(end ,:)';
            SAT_B = 1/2.*H_y_inv_cell{level}*(E_N_y*(u__N-w__0));
        elseif Data.level{bottom_index(1)} > level % if bottom neighbour has higher level (two neighbours)
            Z_bottom_1 = reshape(u((bottom_index(1)-1)*vec_range+1:bottom_index(1)*vec_range),Bs,Bs);
            w__0_1 = Z_bottom_1(end ,:)';
            Z_bottom_2 = reshape(u((bottom_index(2)-1)*vec_range+1:bottom_index(2)*vec_range),Bs,Bs);
            w__0_2 = Z_bottom_2(end ,:)';
            if sum(bottom_neighbours(:,:,1))<sum(bottom_neighbours(:,:,2))
                w__0 = [I_F2C*w__0_1;I_F2C*w__0_2];
            else
                w__0 = [I_F2C*w__0_2;I_F2C*w__0_1];
            end
            SAT_B = 1/2.*H_y_inv_cell{level}*(E_N_y*(u__N-I_single2double'*H_two_blocks*w__0));
        else % if bottom neighbour has lower level
            Z_bottom = reshape(u((bottom_index(1)-1)*vec_range+1:bottom_index(1)*vec_range),Bs,Bs);
            w__0 = Z_bottom(end ,:)';
            w__0_two_blocks = I_single2double*w__0;
            if mod(local_treecode(level),2) == 0 % if this block is left part
                SAT_B = 1/2.*H_y_inv_cell{level}*(E_N_y*(u__N-I_C2F*w__0_two_blocks(1:(Bs+1)/2)));
            else % if this block is right part
                SAT_B = 1/2.*H_y_inv_cell{level}*(E_N_y*(u__N-I_C2F*w__0_two_blocks((Bs+1)/2+1:end)));
            end
            
        end

    
        %u_block_new = (-D_x_cell{level}*u_block+SAT_R+SAT_L);
        u_block_new = (-D_x_cell{level}*u_block-D_y_cell{level}*u_block+SAT_R+SAT_L-SAT_T-SAT_B);
        %u_block_new = (D_y_cell{level}*u_block-SAT_T-SAT_B);
        u_new((block-1)*vec_range+1:block*vec_range) = u_block_new;
    end
    

    
end
    
    
    

