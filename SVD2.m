% ALGORITMO SVD DE GOLUB-KAHAN (2)

function sigma = SVD2(A)

    % reducir A a forma bidiagonal
    J = bidiag(A)
    %J = abs(J);

    % calcular matriz K
    K = J'*J;

    % calcular autovalores de K
    lambda = QRrayleigh(K);

    % vvss de A = raíz cuadrada de autovalores de K
    sigma = sort(sqrt(lambda), "descend");

end