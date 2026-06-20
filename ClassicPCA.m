% ALGORITMO PCA

function [loadings, scores, PVE, PVE_total] = ClassicPCA(X)
    [m,n] = size(X); 
    
    % suponemos X matriz de datos (obs x var) centrada y estandarizada

    % SVD truncada (ahorra coste computacional)
    [U,S,V] = svd(X, 'econ');

    % puntuaciones de cada observación (fila): filas de U*S
    scores = U*S;

    % vectores de cargas de las CP: columnas de V 
    loadings = V;

    % varianza explicada por cada CP: autovalor de X'X = sigma^2/(m-1)
    var_expl = (diag(S).^2)/(m-1);
    
    % porcentaje de varianza explicada
    PVE = (var_expl/sum(var_expl))*100;
    PVE_total = cumsum(PVE);

end

% se ha comprobado que, con una matriz de datos centrada y estandarizada,
% los resultados de esta implementación de PCA son idénticos a los de la
% función pca() de MATLAB

