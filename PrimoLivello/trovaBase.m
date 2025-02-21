function [B] = trovaBase(A)
%TROVABASE Trova una base, ovvero un insieme di indici B tali per cui
%A(B,:) è una matrice invertibile.
%   Data una matrice A, che rappresenta un insieme di vincoli per un
%   problema di programmazione lineare, trova una base, ovvero un
%   sottoinsieme di righe di A, sia esso A(B,:), tale per cui la
%   matrice formata da queste righe è invertibile. Per essere vero ciò,
%   intanto A deve avere un numero di righe maggiore o uguale delle
%   colonne, altrimenti non esisterebbe un sottoinsieme delle righe che
%   forma una matrice quadrata, e poi deve esistere un sottoinsieme di
%   n righe che sono indipendenti tra loro, dove n è il numero di
%   colonne di A. Ciò equivale a dire che il rango di A deve essere
%   pari al numero di colonne. Se entrambe queste affermazioni sono
%   vere, allora esiste almeno una base. L'algoritmo cicla tra tutte le
%   basi, ovvero tutte le possibili combinazioni, in ordine finché non
%   si trova una base la cui matrice di base è invertibile. Siamo
%   sicuri che l'algoritmo si fermerà prima di arrivare all'ultima 
%   base, ovvero [ ... n-2 n-1 n] perché sappiamo che il rango di A è
%   uguale ad n, quindi non c'è bisogno di fare controlli di stop. Per
%   questo non si segnalano errori alla fine del for loop.
%   PARAMETRI
%   A: matrice dei vincoli
%   OUTPUT
%   B: una base per un problema di programmazione lineare con A matrice
%   dei vincoli
%   ESEMPIO
%   A = [0 0.6 0.8;
%       -1 2 0;
%       1 0 -1;
%       -1 0 0;
%       0 -1 0;
%       0 0 -1];
%   trovaBase(A);
if ~exist("A","var")
    error("Errore. La matrice dei vincoli A deve essere passata come parametro.")
end
controlloA(A);
[nVincoli, nVariabili] = size(A);
if nVincoli < nVariabili
    error("Errore. La matrice A deve avere un numero di righe maggiore o uguale al numero delle colonne.");
end
A = sym(A);
if rank(A) ~= nVariabili
    error("Errore. La matrice A non ha rango massimo, quindi non esistono basi.")
end

% Inizializza la combinazione iniziale (1, 2, 3, ..., nVariabili)
B = 1:nVariabili;

while true
    if rank(A(B, :)) == nVariabili
        % Trovata una base valida
        return;
    end
    
    % Genera la prossima combinazione in ordine lessicografico:
    i = nVariabili;
    while i > 0 && B(i) == i - 1 + nVincoli - nVariabili + 1
        i = i - 1;
    end
    if i == 0
        error("Non è stata trovata alcuna base valida.");
    end
    B(i) = B(i) + 1;
    for j = i+1:nVariabili
        B(j) = B(j-1) + 1;
    end
end
end

