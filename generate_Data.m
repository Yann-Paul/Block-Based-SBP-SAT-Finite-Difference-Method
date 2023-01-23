function Data = generate_Data(Domain,Bs,treecode)
    Data = struct('treecode',[],'level',[],'coords_origin',[],'coords_spacing',[],'X',[],'Y',[],'Z',[],'K',[]);
    Data.treecode = treecode;
    
    x_len = Domain(1,2) - Domain(1,1);
    y_len = Domain(2,2) - Domain(2,1);
    
    n_blocks = size(Data.treecode);
    n_blocks = n_blocks(2);
    
    
    for i=1:n_blocks
        Data.level{i} = length(Data.treecode{i});
        treecode = Data.treecode{i};
        coords_origin = [0,0];
        for j=1:length(treecode)
            coords_origin(1) = coords_origin(1) + mod(treecode(j),2)/2^j*x_len;
            coords_origin(2) = coords_origin(2) + (treecode(j)>1)/2^j*y_len;
        end
        Data.coords_origin{i} = coords_origin;
        Data.coords_spacing{i} = [x_len/((Bs-1)*2^(length(treecode))),y_len/((Bs-1)*2^(length(treecode)))];
        x = linspace(Data.coords_origin{i}(1),Data.coords_origin{i}(1)+(Bs-1)*Data.coords_spacing{i}(1),Bs);
        %y = linspace(Data.coords_origin{i}(2),Data.coords_origin{i}(2)-(Bs-1)*Data.coords_spacing{i}(2),Bs);
        y = linspace(Data.coords_origin{i}(2),Data.coords_origin{i}(2)+(Bs-1)*Data.coords_spacing{i}(2),Bs);
        [Data.X{i},Data.Y{i}] = meshgrid(x,y);
        Data.K{i} = zeros(Bs,Bs,5);
    end
    
end