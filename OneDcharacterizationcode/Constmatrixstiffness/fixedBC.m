function [K,f,vec_num] = fixedBC(K,f,vec_num,node_num,elem_num)

    K = [K; zeros(1,vec_num)];
    K = [K, zeros(vec_num + 1,1)];
    vec_num = vec_num + 1;
    K(1 , vec_num) = 1;
    K(vec_num , :) = K(: , vec_num);
    f = [f,0];

end