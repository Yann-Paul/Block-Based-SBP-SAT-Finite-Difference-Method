function [ neighbours_top_array, neighbours_bottom_array,neighbours_top_index,neighbours_bottom_index] = create_bottom_neighbours_list(treecode_array,Data)
    nblocks = size(Data.treecode,2);
    maxlevel = max(cell2mat(Data.level));

    neighbours_bottom_array = -ones(nblocks,maxlevel,2);
    neighbours_top_array = -ones(nblocks,maxlevel,2);
    neighbours_bottom_index = zeros(nblocks,2);    
    neighbours_top_index = zeros(nblocks,2);

    
    % create list of bottom neighbours
    for k = 1:nblocks
        local_treecode = treecode_array(k,:);
        level = Data.level{k};
    

        local_tree_new = transform_bottom(local_treecode,level);
        ref_treecodes = [];
        ref_treecodes(1,:,:) = (local_tree_new.*ones(maxlevel,4)')';

        ref_treecodes(1,level,2) = -1;
    
        if level ~= maxlevel
            ref_treecodes(1,level+1,3) = 2;
            ref_treecodes(1,level+1,4) = 3;

        else 
            ref_treecodes(1,:,3) = -1;
            ref_treecodes(1,:,4) = -1;
        end

        index_neighbours_bottom = sum(prod(ref_treecodes == treecode_array,2),3);
        index_neighbours_bottom = logical(index_neighbours_bottom');
        
        block_neighbours_bottom = treecode_array( index_neighbours_bottom,:);
        index = find( index_neighbours_bottom);
    
        for neighbour = 1:size(block_neighbours_bottom,1)
            neighbours_bottom_array(k,:,neighbour) = block_neighbours_bottom(neighbour,:);
            neighbours_bottom_index(k,neighbour) = index(neighbour);
        end
    end
    
    % create array of top neighbours
    for k = 1:nblocks
        local_treecode = treecode_array(k,:);
        index_neighbours_top = logical(sum(prod(neighbours_bottom_array == local_treecode,2),3));
        block_neighbours_top = treecode_array(index_neighbours_top,:);
        index = find(index_neighbours_top);
    
        for neighbour = 1:size(block_neighbours_top,1)
            neighbours_top_array(k,:,neighbour) = block_neighbours_top(neighbour,:);
            neighbours_top_index(k,neighbour) = index(neighbour);
        end
    end

end