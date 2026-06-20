% HOJA 2, EJERCICIO 2

function [Q,R] = QRhess(A)
    A = hessenberg(A);
    n = size(A, 1);
    R = A;
    Q = eye(n); % inicializamos Q como I
    for i = 1:n-1
        % reflexiones de Householder
        v = R(i:i+1, i);
        if v(1) == 0
            v(1) = v(1) + norm(v);
        else
            v(1) = v(1) + sign(v(1)) * norm(v);
        end
        R(i:i+1, i:n) = R(i:i+1, i:n) - (2 / (v' * v)) * (v * (v' * R(i:i+1, i:n)));
        Q(i:i+1,:) = Q(i:i+1,:) - (2 / (v' * v) * (v * ( v' * Q(i:i+1,:))));
    end    
    Q = Q';