% ALGORITMO ThreSPCA para calcular Sparse PCA
% input:
% X: matriz de datos mxn
% k: parámetro de dispersión (numero de pcs no nulas)
% e: parámetro de error
% ncp: numero de componentes principales dispersas a calcular
% output:
% sparse_loadings: vectores de carga de las sparse CP
% scores: puntuaciones de las observaciones
% PVE: porcentaje de varianza explicada por cada CP
% nnz_z: numero de elementos no nulos de cada CP

function [Z, nnz_z, scores, PVE_total] = ThreSPCA(X,k,epsilon,ncp)
    X_ini = X;
    [m,n] = size(X); 
    Z = zeros(n,ncp);
    l = ceil(1/epsilon);
    nnz_z = zeros(1, ncp);

    for i = 1:ncp

        % calcular matriz de varianzas covarianzas
        A = (X'*X)/(m-1);

        % calcular primeros l vectores y valores singulares de A
        [Ul, Sl,~] = svds(A,l);

        % calcular normas al cuadrado de filas de Ul
        r = sum(Ul.^2, 2);

        % seleccionar filas con norma por encima del umbral y construir
        % matriz
        R = find(r >= epsilon^2 / k);
        SR = eye(n);
        SR = SR(:,R);

        % calcular primer vector singular derecho
        [~, ~, y] = svds(Sl.^(1/2) * Ul' * SR, 1);

        % definir primera CP
        z = SR*y;
        nnz_z(i) = nnz(z);

        % iterar
        [~,~,v] = svds(X,1);
        X = X - X*(v*v');
        Z(:,i) = z;

    end

    scores = X_ini*Z*pinv(Z'*Z);
    TP = scores*Z';
    PVE_total = 100*(trace(TP'*TP)/trace(X_ini'*X_ini));

end


