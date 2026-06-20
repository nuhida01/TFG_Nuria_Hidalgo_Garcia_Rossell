% ALGORITMO QR CON TRASLACON DE RAYLEIGH
% para matriz tridiagonal simétrica

function lambda = QRrayleigh(A)
    n = size(A, 1);           
    lambda = zeros(n, 1); % almacenará los autovalores aproximados de A
    
    if n == 1 % caso base
        lambda(1) = A(1, 1);
    
    elseif n == 2 % caso base
        tr = A(1,1) + A(2,2);
        det = A(1,1)*A(2,2) - A(1,2)*A(2,1);
        lambda(1) = (tr + sqrt(tr^2 - 4*det)) / 2;
        lambda(2) = (tr - sqrt(tr^2 - 4*det)) / 2;
    
    else % algoritmo recursivo
        l_diag = abs(diag(A(1:n, 1:n), -1));
        l_diag_zeros = find(l_diag < eps);
        
        if ~isempty(l_diag_zeros) % matriz reducible
            l_diag_zeros(end+1) = n;
            m = 1;
            
            for i = l_diag_zeros % división por cajas
                sub_A = A(m:i, m:i);
                lambda(m:i) = QRrayleigh(sub_A); % recursión
                m = i + 1;
            end
        
        else % matriz irreducible
            % traslación de Rayleigh
            sk = A(n, n);
            [Q, R] = QRhess(A - sk * eye(n));
            A = R * Q + sk * eye(n);
            lambda = QRrayleigh(A); % recursión
        
        end
    
    end

end