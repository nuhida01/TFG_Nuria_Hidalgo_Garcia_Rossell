% reduccion a forma bidiagonal con matrices de householder

function B = bidiag(A)
    m = size(A, 1);
    n = size(A, 2);

    % inicializamos B como A
    B = A;

    for k = 1:n
        % tomar elementos subdiagonales de la columna k
        x = B(k:m, k);
        % construir la reflexion de householder
        v = x;
        v(1) = v(1) + sign(x(1)) * norm(x);
        v = v / norm(v);
        % aplicar la reflexion de householder por la izquierda
        B(k:m, k:n) = B(k:m, k:n) - 2 * v * (v' * B(k:m, k:n));
        
        if (k < n)
            % tomar elementos de la fila k a la derecha de la diagonal
            x = B(k,k+1:n);
            % construir la reflexion de householder
            v = x';
            v(1) = v(1) + sign(x(1)) * norm(x);
            v = v / norm(v);
            % aplicar la reflexion de householder por la derecha
            B(k:m, k+1:n) = B(k:m, k+1:n) - 2 * (B(k:m, k+1:n) * v) * v';
        end

    end

end