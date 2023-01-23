function [neighbours_right_array,neighbours_left_array,neighbours_right_index,neighbours_left_index] = create_neighbours_list(treecode_array,Data)
    nblocks = size(Data.treecode,2);
    maxlevel = max(cell2mat(Data.level));

    neighbours_right_array = -ones(nblocks,maxlevel,2);
    neighbours_left_array = -ones(nblocks,maxlevel,2);
    neighbours_right_index = zeros(nblocks,2);
    neighbours_left_index = zeros(nblocks,2);
    
    % create list of right neighbours
    for k = 1:nblocks
        local_treecode = treecode_array(k,:);
        level = Data.level{k};
    

        new_tree= transform_right(local_treecode,level);
        ref_treecodes = [];
        ref_treecodes(1,:,:) = (new_tree.*ones(maxlevel,4)')';
        ref_treecodes(1,level,2) = -1;

        if level ~= maxlevel
            ref_treecodes(1,level+1,3) = 0;
            ref_treecodes(1,level+1,4) = 2;
        else
            ref_treecodes(1,:,3) = -1;
            ref_treecodes(1,:,4) = -1;
        end

        index_neighbours_right = sum(prod(ref_treecodes == treecode_array,2),3);
        index_neighbours_right = logical(index_neighbours_right');

        block_neighbours_right = treecode_array(index_neighbours_right,:);
        index = find(index_neighbours_right);
    
        for neighbour = 1:size(block_neighbours_right,1)
            neighbours_right_array(k,:,neighbour) = block_neighbours_right(neighbour,:);
            neighbours_right_index(k,neighbour) = index(neighbour);
        end
    end

    % create array of left neighbours 
    
    for k = 1:nblocks
        local_treecode = treecode_array(k,:);
        index_neighbours_left = logical(sum(prod(neighbours_right_array == local_treecode,2),3));
        block_neighbours_left = treecode_array(index_neighbours_left,:);
        index = find(index_neighbours_left);
    
        for neighbour = 1:size(block_neighbours_left,1)
            neighbours_left_array(k,:,neighbour) = block_neighbours_left(neighbour,:);
            neighbours_left_index(k,neighbour) = index(neighbour);
        end
    end
end