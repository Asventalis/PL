function [] = controlloAc(A,c)
    % CONTROLLOAC Lancia un errore se il vettore c non può essere un vettore
    % che rappresenta la funzione obiettivo per un problema di
    % programmazione lineare in formato primale standard con matrice dei
    % vincoli A.
    %   Controlla le seguenti cose:
    %    - A è passato come argomento ed è una matrice dei vincoli
    %    - c non è vuoto
    %    - c è numerico o simbolico
    %    - c è un vettore colonna
    %    - A ha un numero di colonne pari a length(c)
    %   PARAMETRI
    %   A: matrice dei vincoli
    %   c: vettore della funzione obiettivo
    
    if ~exist("A","var")
        error("Errore. La matrice dei vincoli A va passata come primo argomento.");
    end
    controlloA(A);
    if ~exist("c", "var")
        error("Errore. Il vettore della funzione obiettivo va passato come secondo argomento.");
    end
    if isempty(c)
        error("Errore. Il vettore c è vuoto.");
    end
    controlloNumSym(c);
    if ~iscolumn(c)
        error("Errore. c = %s deve essere un vettore colonna.", stampaInLinea(c));
    end
    [~, nVariabili] = size(A);
    if nVariabili ~= length(c)
        error("Errore. Il numero di variabili (colonne di A) dev'essere uguale alla lunghezza di c.");
    end
end

