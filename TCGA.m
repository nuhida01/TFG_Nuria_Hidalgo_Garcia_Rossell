% IDENTIFICACIÓN DE FIRMA GÉNICA EN CÁNCER DE MAMA

%% CARGA DE DATOS

% datos de expresión génica

data = readtable('brca_tcga_pan_can_atlas_2018\data_mrna_seq_v2_rsem.txt',... 
 'Delimiter', '\t');
head(data,5)

% eliminamos genes sin etiquetar
ids_empty = string(data{:,1}) == '';
data = data(~ids_empty,:);
X = table2array(data(:,3:end));
genes = table2cell(data(:,1));

% trasponemos la matriz para tener filas = muestras, columnas = genes
X = X';
[m,n] = size(X);
disp("Tamaño original de la matriz de datos (muestras x genes): ")
disp(size(X))

% m filas: muestras (pacientes - tumores)
% n columnas: genes

% datos clinicos

clinical_data = readtable('brca_tcga_pan_can_atlas_2018\data_clinical_patient.txt', ...
    'FileType', 'text', ...
    'Delimiter', '\t', ...
    'CommentStyle', '#');

% ver columnas
disp(clinical_data.Properties.VariableNames)

% buscar columna de subtipos de tumores
clinical_data(:, {'PATIENT_ID', 'SUBTYPE'})

% mapeo muestras - subtipos
sample_ids = data.Properties.VariableNames(3:end)';

% normalizar formatos (coger primeros 12 caracteres, puntos a guiones)
patient_ids = cellfun( @(s) strrep(s(1:12), '_', '-'), sample_ids, ...
    'UniformOutput', false);

% buscar en clinical_data
[~, loc] = ismember(patient_ids, clinical_data.PATIENT_ID);

% guardar subtipo de cada muestra
subtypes = repmat({'Unknown'}, length(sample_ids), 1);
found = loc > 0;
subtypes(found) = clinical_data.SUBTYPE(loc(found));

disp('Muestras con subtipo asignado:');
disp(sum(found))


%% LIMPIEZA Y PREPROCESADO DE DATOS

% guardamos los indices para luego saber a que genes corresponden
indices_conservados = 1:size(X, 2);

% tratamiento de datos faltantes
cols_nan = any(isnan(X), 1);
n_cols_nan = sum(cols_nan);
disp('Columnas NaN: ')
disp(n_cols_nan)
X = X(:,~cols_nan);
indices_conservados = indices_conservados(~cols_nan);

rows_nan = any(isnan(X), 2);
n_rows_nan = sum(rows_nan);
disp('Filas NaN: ')
disp(n_rows_nan)
X = X(~rows_nan,:);
% actualizamos la tabla de subtipos de las muestras
subtypes = subtypes(~rows_nan);

% eliminamos genes poco expresados
pct_ceros_cols = sum(X==0 , 1)/m;
cols_ceros = pct_ceros_cols > 0.7;
num_cols_cero = sum(cols_ceros);
disp("número de columnas a eliminar: ")
disp(num_cols_cero)
indices_conservados = indices_conservados(~cols_ceros);

X = X(:,~cols_ceros);
[m,n] = size(X);
disp("tamaño final de la matriz de datos: ")
disp(size(X))

%% CENTRADO Y ESTANDARIZADO DE LA MATRIZ DE DATOS

% centrar la matriz de datos
mu = mean(X,1); % media de cada columna (variable)
X = X - mu; % matriz centrada

% estandarizar
std_genes = std(X, 0, 1);
std_genes(std_genes == 0) = 1;
X_std = X ./ std_genes;

disp("X centrada")


%% PCA
[loadings, scores, PVE, PVE_total] = ClassicPCA(X_std);
disp("Porcentaje de varianza explicada por cada CP")
display(PVE(1:10))
disp("Porcentaje de varianza acumulada para cada CP")
display(PVE_total(1:10))
display(PVE_total(100))
display(PVE_total(200))

% BIPLOT
scores_biplot = scores(:,1:3);
loadings_biplot = loadings(:,1:3);

% colores segun subtipo PAM50 
subtypes_list = {'BRCA_LumA','BRCA_LumB','BRCA_Her2','BRCA_Basal',...
    'BRCA_Normal','Unknown'};
labels_display = {'Luminal A', 'Luminal B', 'HER2', 'Basal', ...
    'Normal', 'Unknown'};
colores = [0.15 0.45 0.90;
           0.10 0.75 0.30;
           1.00 0.55 0.00;
           0.85 0.10 0.10;
           0.60 0.60 0.60;
           0.85 0.85 0.85];

color_muestras = zeros(length(subtypes), 3);
for i = 1:length(subtypes_list)
    idx = strcmp(subtypes, subtypes_list{i});
    color_muestras(idx,:) = repmat(colores(i,:), sum(idx), 1);
end

%  BIPLOT 2D 
figure('Name', 'PCA Biplot 2D');
hold on;

% scores coloreados por subtipo
h_leg = gobjects(length(subtypes_list), 1);
for i = 1:length(subtypes_list)
    idx = strcmp(subtypes, subtypes_list{i});
    h_leg(i) = scatter(scores_biplot(idx,1), scores_biplot(idx,2), ...
                       20, colores(i,:), 'filled', 'MarkerFaceAlpha', 0.6);
end

xlabel(sprintf('CP1 (%.1f%%)', PVE(1)));
ylabel(sprintf('CP2 (%.1f%%)', PVE(2)));
title('PCA: Proyección de las observaciones sobre CP1 y CP2');
legend(h_leg, labels_display, 'Location', 'best');
grid on;
hold off;

%  BIPLOT 3D 
figure('Name', 'PCA Biplot 3D');
hold on;

% scores coloreados por subtipo
h_leg3 = gobjects(length(subtypes_list), 1);
for i = 1:length(subtypes_list)
    idx = strcmp(subtypes, subtypes_list{i});
    h_leg3(i) = scatter3(scores_biplot(idx,1), scores_biplot(idx,2), ...
        scores_biplot(idx,3), 20, colores(i,:), 'filled', ...
        'MarkerFaceAlpha', 0.6);
end

xlabel(sprintf('CP1 (%.1f%%)', PVE(1)));
ylabel(sprintf('CP2 (%.1f%%)', PVE(2)));
zlabel(sprintf('CP3 (%.1f%%)', PVE(3)));
title('PCA: Proyección de las observaciones sobre CP1, CP2 y CP3');
legend(h_leg3, labels_display, 'Location', 'best');
grid on; 
view(45, 25);
hold off;


%% REPETIMOS SELECCIONANDO LOS GENES EN FIRMA PAM50
pam50_genes_hugo = {'UBE2T', 'BIRC5', 'NUF2', 'CDC6', 'CCNB1', 'TYMS', ...
                    'MYBL2', 'CEP55', 'MELK', 'NDC80', 'RRM2', 'UBE2C', ...
                    'CENPF', 'PTTG1', 'EXO1', 'ORC6L', 'ANLN', 'CCNE1', ...
                    'CDC20', 'MKI67', 'KIF2C', 'ACTR3B', 'MYC', 'EGFR', ...
                    'KRT5', 'PHGDH', 'CDH3', 'MIA', 'KRT17', 'FOXC1', ...
                    'SFRP1', 'KRT14', 'ESR1', 'SLC39A6', 'BAG1', 'MAPT', ...
                    'PGR', 'CXXC5', 'MLPH', 'BCL2', 'MDM2', 'NAT1', ...
                    'FOXA1', 'BLVRA', 'MMP11', 'GPR160', 'FGFR4', ...
                    'GRB7', 'TMEM45B', 'ERBB2'};

rows_pam50 = ismember(data{:,1}, pam50_genes_hugo);
data_pam50 = data(rows_pam50, :);

X_PAM50 = table2array(data_pam50(:,3:end));
genes_pam50 = data_pam50(:,1);

% trasponemos la matriz para tener filas = muestras, columnas = genes
X_PAM50 = X_PAM50';
[m,n] = size(X_PAM50);
disp("Tamaño original de la matriz de datos PAM50 (muestras x genes): ")
disp(size(X_PAM50))

% m filas: muestras (pacientes - tumores)
% n columnas: genes


%% LIMPIEZA Y PREPROCESADO DE DATOS

% tratamiento de datos faltantes
cols_nan = any(isnan(X_PAM50), 1);
n_cols_nan = sum(cols_nan);
disp('Columnas NaN: ')
disp(n_cols_nan)
X_PAM50 = X_PAM50(:,~cols_nan);

rows_nan = any(isnan(X_PAM50), 2);
n_rows_nan = sum(rows_nan);
disp('Filas NaN: ')
disp(n_rows_nan)
X_PAM50 = X_PAM50(~rows_nan,:);

% eliminación de genes poco expresados
pct_ceros_cols = sum(X_PAM50==0, 1)/m;
cols_ceros = pct_ceros_cols > 0.7;
num_cols_cero = sum(cols_ceros);
disp("número de columnas a eliminar: ")
disp(num_cols_cero)

X_PAM50 = X_PAM50(:,~cols_ceros);
[m,n] = size(X_PAM50);
disp("tamaño final de la matriz de datos: ")
disp(size(X_PAM50))

%% CENTRADO Y ESTANDARIZADO DE LA MATRIZ DE DATOS

% centrar la matriz de datos
mu = mean(X_PAM50,1); % media de cada columna (variable)
X_PAM50 = X_PAM50 - mu; % matriz centrada

% estandarizar
std_genes = std(X_PAM50, 0, 1);
std_genes(std_genes == 0) = 1;
X_PAM50 = X_PAM50 ./ std_genes;

disp("X_PAM50 centrada y estandarizada")


%% PCA
[loadings, scores, PVE, PVE_total] = ClassicPCA(X_PAM50);

disp("Porcentaje de varianza explicada por cada CP")
display(PVE(1:10))
disp("Porcentaje de varianza acumulada para cada CP")
display(PVE_total(1:10))




% BIPLOT
scores_biplot   = scores(:,1:3);
loadings_biplot = loadings(:,1:3);

% escalado para ver las flechas de loadings
scale = max(abs(scores_biplot),[],1) ./ max(abs(loadings_biplot),[],1)*0.7;
loadings_scaled = loadings_biplot .* scale;

color_muestras = zeros(length(subtypes), 3);
for i = 1:length(subtypes_list)
    idx = strcmp(subtypes, subtypes_list{i});
    color_muestras(idx,:) = repmat(colores(i,:), sum(idx), 1);
end

%  BIPLOT 2D 
figure('Name', 'PCA Biplot 2D');
hold on;

% scores coloreados por subtipo
h_leg = gobjects(length(subtypes_list), 1);
for i = 1:length(subtypes_list)
    idx = strcmp(subtypes, subtypes_list{i});
    h_leg(i) = scatter(scores_biplot(idx,1), scores_biplot(idx,2), ...
                       20, colores(i,:), 'filled', 'MarkerFaceAlpha', 0.6);
end

% loadings como flechas (quiver)
n_load = size(loadings_scaled, 1);
quiver(zeros(n_load,1), zeros(n_load,1), ...
       loadings_scaled(:,1), loadings_scaled(:,2), ...
       0, 'k', 'LineWidth', 0.8, 'MaxHeadSize', 0.5);

xlabel(sprintf('PC1 (%.1f%%)', PVE(1)));
ylabel(sprintf('PC2 (%.1f%%)', PVE(2)));
title('PAM50 PCA: Biplot con CP1 y CP2');
legend(h_leg, labels_display, 'Location', 'best');
grid on;
hold off;

% BIPLOT 3D 
figure('Name', 'PCA Biplot 3D');
hold on;

% scores coloreados por subtipo
h_leg3 = gobjects(length(subtypes_list), 1);
for i = 1:length(subtypes_list)
    idx = strcmp(subtypes, subtypes_list{i});
    h_leg3(i) = scatter3(scores_biplot(idx,1), scores_biplot(idx,2), ...
        scores_biplot(idx,3), ...
                         20, colores(i,:), 'filled', 'MarkerFaceAlpha', 0.6);
end

% loadings como flechas (quiver3)
quiver3(zeros(n_load,1), zeros(n_load,1), zeros(n_load,1), ...
        loadings_scaled(:,1), loadings_scaled(:,2), loadings_scaled(:,3), ...
        0, 'k', 'LineWidth', 0.8, 'MaxHeadSize', 0.5);

xlabel(sprintf('PC1 (%.1f%%)', PVE(1)));
ylabel(sprintf('PC2 (%.1f%%)', PVE(2)));
zlabel(sprintf('PC3 (%.1f%%)', PVE(3)));
title('PAM50 PCA: Biplot con CP1, CP2 y CP3');
legend(h_leg3, labels_display, 'Location', 'best');
grid on; 
view(45, 25);
hold off;




%% ThreSPCA
% con la matriz de datos completa X
% solo centrada, no estandarizada para que el umbral funcione
[sparseZ, nnz, sparse_scores, sparse_PVE_total] = ThreSPCA(X,100,0.45,3);

%% PLOT: SCATTER DE SCORES
figure('Name', 'sPCA Scatter Scores 2D');
hold on;
h_leg = gobjects(length(subtypes_list), 1);
for i = 1:length(subtypes_list)
    idx = strcmp(subtypes, subtypes_list{i});
    h_leg(i) = scatter(sparse_scores(idx,1), sparse_scores(idx,2), 20, ...
        colores(i,:), 'filled', 'MarkerFaceAlpha', 0.6);
end
xlabel('Sparse PC1'); ylabel('Sparse PC2');
title('ThreSPCA: Proyección de las observaciones sobre CP1 y CP2');
legend(h_leg, labels_display, 'Location', 'best');
grid on; 
hold off;

%% PLOT: SCATTER 3D DE SCORES
figure('Name', 'sPCA Scatter Scores 3D');
hold on;
h_leg3 = gobjects(length(subtypes_list), 1);
for i = 1:length(subtypes_list)
    idx = strcmp(subtypes, subtypes_list{i});
    h_leg3(i) = scatter3(sparse_scores(idx,1), sparse_scores(idx,2), ...
        sparse_scores(idx,3),20, colores(i,:), 'filled', ...
        'MarkerFaceAlpha', 0.6);
end
xlabel('Sparse PC1'); ylabel('Sparse PC2'); zlabel('Sparse PC3');
title('ThreSPCA: Proyección de las observaciones sobre CP1, CP2 y CP3');
legend(h_leg3, labels_display, 'Location', 'best');
grid on; 
view(45, 25); 
hold off;


%% Intersección con genes de PAM50 para cada Sparse PC
fprintf('\n GENES PAM50 EN CADA SPARSE PC \n');
gene_names_X = genes(indices_conservados);

for pc = 1:size(sparseZ, 2)

    % índices con loading != 0
    nz_idx = find(sparseZ(:, pc) ~= 0);
    nz_genes{pc} = gene_names_X(nz_idx);
    nz_vals = sparseZ(nz_idx, pc);

    % intersección con PAM50
    [genes_en_pam50, ids, ~] = intersect(nz_genes{pc}, pam50_genes_hugo);
    loadings_pam50 = nz_vals(ids);

    fprintf(' Sparse PC%d ', pc);
    fprintf('Genes no nulos totales : %d\n', length(nz_idx));
    fprintf('Genes en PAM50: %d / %d\n', length(genes_en_pam50), length(nz_idx));
    fprintf('Cobertura PAM50: %.1f%%\n\n', length(genes_en_pam50)/50*100);

    % tabla ordenada por abs de loading descendente
    [~, ord] = sort(abs(loadings_pam50), 'descend');
    fprintf('  %-12s  %8s\n', 'Gen', 'Loading');
    fprintf('  %-12s  %8s\n', '------------', '--------');
    for j = 1:length(ord)
        fprintf('  %-12s  %+8.4f\n', genes_en_pam50{ord(j)}, ...
            loadings_pam50(ord(j)));
    end
end
