        function dX =dynamics(X,T,J_inv,J)
            dX = zeros(12,1);
           % o = X(1:2);
           % R_now = reshape(X(1:9), 3, 3);
           % W_now = X(10:12);
            R_now = [1,0,0;0,1,0;0,0,1];
            W_now = [1;0;0];
            dR = R_now*hat_map(W_now);
            dW = J_inv*(-cross(W_now,J*W_now) +T);
            
            dX(1:9) = reshape(dR, 9, 1);
            dX(10:12) = dW;
        end

