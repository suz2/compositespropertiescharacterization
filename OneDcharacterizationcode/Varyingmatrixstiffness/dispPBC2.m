function [K,f,PBC_option] = dispPBC2(K,f,bcnode,vec_num,uload)

    PBC_option = 1;
    
    i = bcnode(1);
    j = bcnode(2);
    (K(i,j) + K(j,i))
    
    for k = 1 : length(K(:,1))
       for l = 1 :  length(K(1,:))
           K(k,l) = K(k,l) - K(k,j) * K(j,l)/(K(i,j) + K(j,i));
       end
       f(k) = f(k) - K(k,j)*(f(j) + K(j,i)*uload)/(K(i,j) + K(j,i));
    end

    K(j,:) = 0;
    K(:,j) = 0;
    K(j,j) = 1;
    f(j) = 0;

end