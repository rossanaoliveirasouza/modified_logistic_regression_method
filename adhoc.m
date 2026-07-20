function [alpha, p] = adhoc(q, A0, A1, k0, k1)
% Define a função "adhoc" que retorna o valor "p".
% Entrada:
%   q  - vetor de consulta, (tipo o paciente que queremos diagnosticar, q de question, pergunta)
%   A0 - conjunto de vetores da classe 0.
%   A1 - conjunto de vetores da classe 1.
%   k0 - número de vizinhos mais próximos a considerar da classe 0.
%   k1 - número de vizinhos mais próximos a considerar da classe 1.
% Saída:
%   p  - probabilidade calculada para o vetor "q".
%   A saída p na função é uma probabilidade logística calculada para o vetor q, que é uma medida da probabilidade de q pertencer à classe 1 (ou à classe positiva).

% Calcula a distância euclidiana entre "q" e todos os vetores em "A0".
dist = pdist2(A0, q); 

% Ordena as distâncias em ordem crescente.
% "val" contém os valores ordenados e "pos" contém os índices dos vetores
% correspondentes em "A0".
[val, pos] = sort(dist);

% Seleciona os "k0" vetores mais próximos de "q" no conjunto "A0" com base nas distâncias calculadas.
M = A0(pos(1:k0),:);

% Calcula a distância euclidiana entre "q" e todos os vetores em "A1".
dist = pdist2(A1, q);

% Ordena as distâncias em ordem crescente para os vetores de "A1".
[val, pos] = sort(dist);

% Adiciona os "k1" vetores mais próximos de "q" no conjunto "A1" ao final da matriz "M", logo após os k0's, que já foram adicionados à M.
M = [M; A1(pos(1:k1), :)];

% Adiciona uma coluna de uns à matriz "M". Isso é necessário para representar o termo constante no modelo linear (adiciona o intercepto).
M = [M ones(k0 + k1, 1)];

% Inicializa o vetor "b" com zeros. Ele terá tamanho "k0 + k1", ou seja, M-1, pois M tem o intercepto.
b = zeros(k0 + k1, 1);

% Define os primeiros "k0" valores de "b" como log(0.001 / (1 - 0.001)). Isso representa um peso fixo para os vizinhos da classe 0.
% Define os últimos "k1" valores de "b" como log(0.999 / (1 - 0.999)). Isso representa um peso fixo para os vizinhos da classe 1.
b(1:k0) = log(0.001/(1 - 0.001));
b(k0+1: end) = log(0.999/(1 - 0.999));

% Calcula o vetor "alpha" resolvendo o sistema linear regularizado: alpha = (I + M'*M) \ (M'*b)
% Aqui, "I" é a matriz identidade e "M'*M" representa a matriz de covariância.  
% alpha = (eye(size(M, 2)) + M'*M)\(M'*b);
alpha = resolve(M,b); % modelo com a matriz esparça
%alpha = M\b;

% Adiciona o termo constante "1" ao vetor "q" e calcula o produto escalar com "alpha". O resultado é armazenado em "aux".
aux = [q 1]*alpha;

% Calcula a exponencial de "aux". Isso é parte da transformação logística.
num = exp(aux);

% Calcula a probabilidade "p" usando a fórmula logística: p = e^(aux) / (1 + e^(aux)). Função inversa da logit.
p = num./(1 + num);
end


  
