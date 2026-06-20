% ALGORITMO SVD DE GOLUB-KAHAN (1)

function sigma = SVD1(A)
    m = size(A,1);
    n = size(A,2);
    
    % reducir A a forma bidiagonal
    J = bidiag(A);
    
    % almacenar valores de la diagonal y superdiagonal de J en vectores
    a = diag(J);
    b = diag(J,1);

    % construir la matriz S tridiagonal especial
    S=zeros(2*n);
    for k = 1:n
        S(2*k, 2*k-1) = a(k);
        S(2*k-1, 2*k) = a(k);
        if k < n
            S(2*k, 2*k+1) = b(k);
            S(2*k+1, 2*k) = b(k);
        end
    end

    % calcular autovalores de S = vvss de A
    lambda = QRrayleigh(S);
    sigma = lambda(lambda>=0);

end
