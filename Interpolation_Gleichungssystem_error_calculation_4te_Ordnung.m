
errors = zeros(27,4);
counter = 1;
%for boundary_height_F2C=2:4;
for boundary_height_F2C=5;

    %for boundary_width_F2C = 13;
    for boundary_width_F2C = 10:13;
        
        %for stencil_width_F2C = 9;
        for stencil_width_F2C = 7:2:11;

            Lap = sym([]);
            equal_con_1 = sym([]);
            equal_con_2 = sym([]);
            equal_con_3 = sym([]);
            equal_con_4 = sym([]);
            equal_con_5 = sym([]);
            equal_con_6 = sym([]);
            equal_con_7 = sym([]);

            if boundary_height_F2C == 5 & boundary_width_F2C <= 11 
                continue
            elseif boundary_height_F2C == 5 & boundary_width_F2C <= 12 & stencil_width_F2C >= 9
                continue
            elseif boundary_height_F2C == 5 & boundary_width_F2C <= 13 & stencil_width_F2C >= 11
                continue
            elseif boundary_height_F2C == 2 & stencil_width_F2C >= 11
                continue
            end

            % if mod(boundary_width_F2C,4) == 3
            %     stencil_height_F2C = boundary_width_F2C+2-2*boundary_height_F2C;
            % elseif mod(boundary_width_F2C,4) == 1
            %     stencil_height_F2C = boundary_width_F2C+4-2*boundary_height_F2C;
            % end
            
            if boundary_width_F2C == 12
                if boundary_height_F2C<=4
                    stencil_height_F2C = 15-2*boundary_height_F2C;
                else
                   stencil_height_F2C = 15-2*4;
                end
            elseif boundary_height_F2C<=4
                stencil_height_F2C = 13-2*boundary_height_F2C;
            else
                stencil_height_F2C = 13-2*4;
            end
            
            boundary_height_C2F = boundary_width_F2C;
            
            rows = boundary_height_F2C*2+stencil_height_F2C;
            midrows = (rows+1)/2;
            cols = 2*rows-1;
            midcols = rows;
            
            I_F2C = sym('A', [rows,cols],'real');
            
            for i=1:boundary_height_F2C
                for j=boundary_width_F2C+1:cols
                    I_F2C(i,j) = 0;
                end
            end
            
            for j=1:(cols-stencil_width_F2C)/2
                I_F2C(midrows,j) = 0;
            end
            
            for j=cols:-1:midcols+1
                I_F2C(midrows,j) = I_F2C(midrows,cols-j+1);
            end
            
            for i=midrows-(stencil_height_F2C-1)/2:midrows-1
                for j=1:cols
                    I_F2C(i,j) = 0;
                end
                for j=i*2-1-(stencil_width_F2C-1)/2:i*2-1+(stencil_width_F2C-1)/2
                %for j=midcols-2*(midrows-i)-(stencil_width_F2C-1)/2:midcols-2*(midrows-i)+(stencil_width_F2C-1)/2
                    I_F2C(i,j) = I_F2C(midrows,j+2*(midrows-i));
                end
            end
            
            
            
            
            
            for i=rows:-1:midrows
                for j=cols:-1:1
                    I_F2C(i,j) = I_F2C(rows-i+1,cols-j+1);
                end
            end
            
            I_F2C
            
            parameters = reshape(I_F2C(1:boundary_height_F2C,1:boundary_width_F2C),1,[]);
            parameters = [parameters,I_F2C(midrows,(cols-stencil_width_F2C)/2+1:midcols)]
            
            % create H-matrix for fine grid 
            h_vec = ones(cols,1);
            h_vec(1) = 17/48;
            h_vec(2) = 59/48;
            h_vec(3) = 43/48;
            h_vec(4) = 49/48;
            h_vec(end:-1:end-3) = h_vec(1:4);
            H_fine = diag(h_vec);
            
            % create H-matrix for coarse grid 
            h_vec = ones(rows,1);
            h_vec(1) = 17/48;
            h_vec(2) = 59/48;
            h_vec(3) = 43/48;
            h_vec(4) = 49/48;
            h_vec(end:-1:end-3) = h_vec(1:4);
            H_coarse = 2*diag(h_vec);
            
            % create coarsening matrix from refinement matrix
            %I_F2C = (H_fine*I_C2F*inv(H_coarse))'
            I_C2F = (H_coarse*I_F2C*inv(H_fine))'
            
            equal_con_0 = sum(I_C2F(1:midcols,:),2) -1;
            equal_con_1 = sum(I_F2C(1:midrows,:),2) -1;
            for i=1:boundary_height_F2C+1
                dist = [1:cols]-(2*i-1);
                equal_con_2(i,1) = sum(I_F2C(i,:).*dist) ;
            end
            for i=1:midcols
                dist = [1:rows]-(i+1)/2;
                equal_con_3(i,1) = sum(I_C2F(i,:).*dist);
            end
            
            % interpolation property for power function x^2 in the interior
            for i=boundary_height_C2F+1:midcols
                dist = ([1:rows]-(i+1)/2).^2;
                equal_con_4(i-boundary_height_C2F,1) = sum(I_C2F(i,:).*dist);
            end
            for i=boundary_height_F2C+1
                dist = ([1:cols]-(2*i-1)).^2;
                equal_con_5(i-boundary_height_F2C,1) = sum(I_F2C(i,:).*dist) ;
            end
            
            % interpolation property for power function x^3 in the interior
            for i=boundary_height_C2F+1:midcols
                dist = ([1:rows]-(i+1)/2).^3;
                equal_con_6(i-boundary_height_C2F,1) = sum(I_C2F(i,:).*dist);
            end
            for i=boundary_height_F2C+1
                dist = ([1:cols]-(2*i-1)).^3;
                equal_con_7(i-boundary_height_F2C,1) = sum(I_F2C(i,:).*dist) ;
            end
            
            constraints = [equal_con_0;equal_con_1;equal_con_2;equal_con_3;equal_con_4;equal_con_5;equal_con_6;equal_con_7];
            
            S = solve(constraints == 0,parameters,'ReturnConditions',true)
            
            for i=1:rows
                for j=1:cols
                    if not(isequal(I_F2C(i,j),sym(0)))
                        I_F2C(i,j) = S.(string(I_F2C(i,j)));
                        %I_F2C(i,j) = S.(tostring(I_F2C(i,j)));
                    end
                end
            end
            I_C2F = (H_coarse*I_F2C*inv(H_fine))'
            
            dist = ([1:rows]'-([1:boundary_height_C2F]+1)/2).^2';
            error1 = sum((I_C2F(1:boundary_height_C2F,:).*dist),2).^2;
            dist = (([1:cols]'-(2.*[1:boundary_height_F2C]-1))/2).^2';
            error2 = sum((I_F2C(1:boundary_height_F2C,:).*dist),2).^2 ;
            error = sum(error1)+sum(error2)
            

            for i=1:length(S.parameters)
                Lap(i) = diff(error,S.parameters(i));
            end
            
            if hasSymType(error,'variable')
                Sol = solve(Lap==0,S.parameters,"Real",true)
            
            
                I_C2F = subs(I_C2F,Sol) 
                I_F2C = subs(I_F2C,Sol)
                error = double(subs(error,Sol));
            end
                
            x_coarse = 0:1/(rows-1):1;
            x_fine = 0:1/(cols-1):1;
            
            err_coarse = double(sum((x_coarse.^2'-I_F2C*x_fine.^2').^2))
            err_fine = double(sum((x_fine.^2'-I_C2F*x_coarse.^2').^2))
            %err_norm = (err_coarse+err_fine)*((rows)^2)
            %err_norm = (err_coarse+err_fine)

            errors(counter,1) = boundary_height_F2C;
            errors(counter,2) = boundary_width_F2C;
            errors(counter,3) = stencil_width_F2C;
            %errors(counter,4) = err_norm;
            errors(counter,4) = double(error);
            counter = counter+1;

        end
    end
end

save("interpolation_errors.mat","errors")

