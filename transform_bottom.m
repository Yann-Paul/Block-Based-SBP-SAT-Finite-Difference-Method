function local_tree_new = transform_bottom(local_tree,lvl)
    local_tree_new = local_tree;
    i = lvl;
    while i>=1 && local_tree_new(i)+2 < 4 
        local_tree_new(i) = local_tree_new(i)+2;
        i = i-1;
    end
    if i>=1
        local_tree_new(i) = local_tree_new(i)-2;
    end
%     for i = lvl:-1:1
%         if local_tree_new(i)+2 < 4
%             local_tree_new(i) = local_tree_new(i)+2;
%         else
%             local_tree_new(i) = local_tree_new(i)-2;
%             break
%         end
%     end
end
