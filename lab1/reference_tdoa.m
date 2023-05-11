function [h]=reference_tdoa(t,x,u,th)
    y_permute = [-1 1 0 0; 
                 -1 0 1 0; 
                 -1 0 0 1];
    h = y_permute*r0_tdoa(0,[x(1),x(2),0],0,th);
end