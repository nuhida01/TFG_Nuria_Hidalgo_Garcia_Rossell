% HOJA 2, EJERCICIO 1

function H = hessenberg(A)
    n = size(A, 1); 
    H = A;
    for i = 1:n-2
        % reflexiones de Householder
        v = H(i+1:n, i);
        if v(1) == 0
            v(1) = v(1) + norm(v);
        else
            v(1) = v(1) + sign(v(1)) * norm(v);
        end        
        H(i+1:n, i:n) = H(i+1:n, i:n) - (2 / (v' * v)) * (v * (v' * H(i+1:n, i:n)));
        H(1:n, i+1:n) = H(1:n, i+1:n) - (2 / (v' * v)) * ((H(1:n, i+1:n) * v) * v');
    end
                          
