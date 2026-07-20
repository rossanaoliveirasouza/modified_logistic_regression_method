% Carregar dados
load("base3/sem_sepse.mat");
load("base3/com_sepse.mat");

% Converter tabelas em arrays
pss = table2array(pacientes_sem_sepse_base3);
pcs = table2array(pacientes_com_sepse_base3);
genes = pacientes_com_sepse_base3.(1);
genes = string(genes);

% Obter tamanhos
[linhas_s, colunas_s] = size(pss);
[linhas_c, colunas_c] = size(pcs);

% Inicializar arrays para armazenar os dados convertidos
Pss = zeros(linhas_s, colunas_s);
Pcs = zeros(linhas_c, colunas_c);

% Converter pss de char para numérico
for i = 1:linhas_s
    for j = 1:colunas_s
        aux = pss{i,j};
        aux = char(aux);
        Pss(i,j) = str2double(aux); % Usar str2double para tratar NaN corretamente
    end
end

% Converter pcs de char para numérico
for i = 1:linhas_c
    for j = 1:colunas_c
        aux = pcs{i,j};
        aux = char(aux);
        Pcs(i,j) = str2double(aux);
    end
end

% Remover a coluna de genes dos arrays, já está mapeado na variável genes
Pss = [Pss(:,2:end)];
Pcs = [Pcs(:,2:end)];

A = [Pss,Pcs];

[T,S,V] = svd(A, 'econ'); % OU [T,S,V] = svd(A);
s = diag(S); % ver quantos grupos tem
s = s*(1/sum(s));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w');

plot(1:length(s), s, '-o', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 5);
hold on;

[~, idx] = max(abs(diff(diff(s))));
idx = idx + 1;

plot(idx, s(idx), 'o', ...
    'MarkerSize', 9, ...
    'LineWidth', 1.8);

text(idx, s(idx), sprintf('  Cotovelo: componente %d', idx), ...
    'FontSize', 10, ...
    'VerticalAlignment', 'bottom');

xlabel('Componente singular');
ylabel('Valor singular normalizado');
title('Espectro singular da matriz de expressão gênica' ,'antes a seleção de atributos (GSE69063)', 'FontSize', 14, 'FontWeight', 'bold');

grid on;
box on;

ax = gca;
ax.FontSize = 11;
ax.GridAlpha = 0.18;

xlim([1 length(s)]);
ylim([0 max(s)*1.10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;

% Define cores
burntOrange = [0.80 0.33 0.00]; % Laranja queimado (Sepse)
greenColor = [0.33 0.80 0.33]; % Verde (Grupo controle)

% Calcula os pontos no espaço 3D
Aux = S * V'; % Combinação dos padrões => A = T * (S * V^T)
x = Aux(1,:);
y = Aux(2,:);
z = Aux(3,:);

% Plota esferas para o grupo controle (1:33) - Verde
scatter3(x(1:33), y(1:33), z(1:33), 100, greenColor, 'filled', ...
         'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 1);

% Plota esferas para os pacientes com sepse (34:90) - Laranja queimado
scatter3(x(34:90), y(34:90), z(34:90), 100, burntOrange, 'filled', ...
         'MarkerEdgeColor', [0.3 0.3 0.3], 'MarkerFaceAlpha', 1);

% Configura a visualização 3D
view(3);
grid on;
lighting gouraud; % Suaviza a iluminação
camlight('right'); % Adiciona luz lateral

% Adiciona título e rótulos
title('NCBI GSE69063 - Singular Value Decomposition (SVD)', ...
    'FontSize', 14, 'FontWeight', 'bold');
%title('Pacients Space of GSE69063' ,'Contribution of Singular Vectors', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('First Singular Vector Coefficient', 'FontSize', 12);
ylabel('Second Singular Vector Coefficient', 'FontSize', 12);
zlabel('Third Singular Vector Coefficient', 'FontSize', 12);

% Adiciona legenda na melhor posição
legend({'Control group', 'Sepsis group'}, ...
       'FontSize', 12, 'Location', 'best');
hold off;


b(1:33) = log(0.00001 / (1 - 0.00001));
b(34:90) = log(0.999999 / (1 - 0.99999));
[alpha,x] = resolve(A',b);
[valores, pos] = sort(alpha);

Ar = [A([pos(1:33) pos(end-32:end)], :)];


[T,S,V] = svd(Ar, 'econ'); % OU [T,S,V] = svd(A);
s = diag(S); % ver quantos grupos tem
s = s*(1/sum(s));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w');

plot(1:length(s), s, '-o', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 5);
hold on;

[~, idx] = max(abs(diff(diff(s))));
idx = idx + 1;

plot(idx, s(idx), 'o', ...
    'MarkerSize', 9, ...
    'LineWidth', 1.8);

text(idx, s(idx), sprintf('  Cotovelo: componente %d', idx), ...
    'FontSize', 10, ...
    'VerticalAlignment', 'bottom');

xlabel('Componente singular');
ylabel('Valor singular normalizado');
title('Espectro singular da matriz de expressão gênica' ,'após a seleção de atributos (GSE69063)', 'FontSize', 14, 'FontWeight', 'bold');

grid on;
box on;

ax = gca;
ax.FontSize = 11;
ax.GridAlpha = 0.18;

xlim([1 length(s)]);
ylim([0 max(s)*1.10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;

% Calcula os pontos no espaço 3D
Aux = S * V'; % Combinação dos padrões => A = T * (S * V^T)
x = Aux(1,:);
y = Aux(2,:);
z = Aux(3,:);

% Plota esferas para o grupo controle (1:33) - Verde
scatter3(x(1:33), y(1:33), z(1:33), 100, greenColor, 'filled', ...
         'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 1);

% Plota esferas para os pacientes com sepse (34:90) - Laranja queimado
scatter3(x(34:90), y(34:90), z(34:90), 100, burntOrange, 'filled', ...
         'MarkerEdgeColor', [0.3 0.3 0.3], 'MarkerFaceAlpha', 1);

% Configura a visualização 3D
view(3);
grid on;
lighting gouraud; % Suaviza a iluminação
camlight('right'); % Adiciona luz lateral

% Adiciona título e rótulos
%title('Patients Space of GSE69063 (After Feature Selection)' ,'Contribution of Singular Vectors', 'FontSize', 14, 'FontWeight', 'bold');
title({'NCBI GSE69063 - Singular Value Decomposition (SVD)' ...
    'After the feature selection'},...
    'FontSize', 14, 'FontWeight', 'bold');
xlabel('First Singular Vector Coefficient', 'FontSize', 12);
ylabel('Second Singular Vector Coefficient', 'FontSize', 12);
zlabel('Third Singular Vector Coefficient', 'FontSize', 12);

% Adiciona legenda na melhor posição
legend({'Control group', 'Sepsis group'}, ...
       'FontSize', 12, 'Location', 'best');
hold off;

clearvars