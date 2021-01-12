function [K,f] = globalstiffness(elem_num,node_num,le,dT,E_matrix,E_fiber,alpha_matrix,alpha_fiber,material_no,vec_num)

    % initialize global stiffness and force vector
    K = zeros(node_num, node_num);
    f = zeros(1, node_num);

    for i = 1: elem_num

        % construct element stiffness matrix Ke and force vector fe
%         Ke = [ 1 -1 ; -1 1 ] / le(i);
%         fe = [-1 1] * dT;

        if material_no(i) == 1
             E = E_fiber;
             alpha = alpha_fiber;
        else
            E = E_matrix(i);
            alpha = alpha_matrix(i);
        end
        
        Ke = [ 1 -1 ; -1 1 ] / le(i) * E;
        fe = [-1 1] * dT * E * alpha;

        % Construct global stiffness matrix and force vector
        if i == 1
            K( i : (i + 1) , i : (i + 1)) = Ke;
            f(i : (i + 1)) = fe;
        elseif i > 1
            K( i : (i + 1) , i : (i + 1)) = K( i : (i + 1) , i : (i + 1) ) + Ke;
            f(i : (i + 1)) = f(i : (i + 1)) + fe;
        end

    end
    
end