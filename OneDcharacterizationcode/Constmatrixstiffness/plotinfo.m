function [l,lm,m,alpha_matrix_plot,E_matrix_plot] = plotinfo(node_num,elem_num,material_no,alpha_matrix_elem,E_matrix_elem,le,l)
    
    % nodal coordinates
    l(1) = 0;
    for i = 2 : node_num
       l(i) = l(i-1) + le(i-1);
    end

    % plotting material assignment
    for i = 1 : elem_num
        
       lm(2*i-1) = l(i);
       lm(2*i) = l(i+1);
       m(2*i-1) = material_no(i);
       m(2*i) = material_no(i);
       alpha_matrix_plot(2*i-1) = alpha_matrix_elem(i);
       alpha_matrix_plot(2*i) = alpha_matrix_elem(i);
       E_matrix_plot(2*i-1) = E_matrix_elem(i);
       E_matrix_plot(2*i) = E_matrix_elem(i);
       
    end

end