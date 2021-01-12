function [K,f,vec_num,PBCoption] = dispPBC(K,f,bcnode,vec_num,uload)

    PBCoption = 1;

    K = [K; zeros(1,vec_num)];
    K = [K, zeros(vec_num + 1,1)];
    vec_num = vec_num + 1;

    K(bcnode(1) , vec_num) = 1;
    K(bcnode(2) , vec_num) = -1;
    K(vec_num , :) = K(: , vec_num);
    f = [f;-uload];

end