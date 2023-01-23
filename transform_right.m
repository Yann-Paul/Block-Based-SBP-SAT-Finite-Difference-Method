function new_tree = transform_right_new(local_tree,lvl) 
    new_tree = local_tree;
    i = lvl;
    while i >= 1 && mod(new_tree(i),2) == 1
        new_tree(i) = new_tree(i) - 1;
        i = i-1;
    end
    if i >= 1
        new_tree(i) = new_tree(i) + 1;
    end

end
