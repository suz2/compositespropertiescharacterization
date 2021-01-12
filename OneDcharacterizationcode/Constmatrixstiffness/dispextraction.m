function [l,u_fiber,l_fiber] = dispextraction(node_num,elem_num,material_no,u,le,l)
    
    % nodal coordinates
    l(1) = 0;
    for i = 2 : node_num
       l(i) = l(i-1) + le(i-1);
    end

    % get mid-fiber displacement
    j = 0;
    for i = 1 : elem_num
        if material_no(i) == 1

            if (i - 1) == 0
                u_fiber_0 = u(i);
                l_fiber_0 = l(i);
            else
                if material_no(i - 1) == 0
                    u_fiber_0 = u(i);
                    l_fiber_0 = l(i);
                end
            end

            if (i + 1) == node_num
                u_fiber_1 = u(i + 1);
                l_fiber_1 = l(i + 1);
                j = j + 1;
                u_fiber(j) = (u_fiber_0 + u_fiber_1)/2;
                l_fiber(j) = (l_fiber_0 + l_fiber_1)/2; 

            else

                if material_no(i + 1) == 0
                    u_fiber_1 = u(i + 1);
                    l_fiber_1 = l(i + 1);

                    j = j + 1;
                    u_fiber(j) = (u_fiber_0 + u_fiber_1)/2;
                    l_fiber(j) = (l_fiber_0 + l_fiber_1)/2; 

                end
            end

        end

    end
    
end