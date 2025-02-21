function [] = controlloAb(A,b)
% CONTROLLOAB Lancia un errore se il vettore b non può essere un vettore
% dei termini noti per un problema di programmazione lineare in formato
% primale standard con matrice dei vincoli A.
%   Controlla le seguenti cose:
%    - A è passato come argomento ed è una matrice dei vincoli
%    - b non è vuoto
%    - b è numerico o simbolico
%    - b è un vettore colonna
%    - A ha un numero di righe pari a length(b)
%   PARAMETRI
%   A: matrice dei vincoli
%   b: vettore dei termini noti
    if ~exist("A","var")
        error("Errore. La matrice dei vincoli A va passata come primo argomento.");
    end
    controlloA(A);
    if ~exist("b", "var")
        error("Errore. Il vettore dei termini noti b va passato come secondo argomento.");
    end
    if isempty(b)
        error("Errore. Il vettore b è vuoto.");
    end
    controlloNumSym(b);
    if ~iscolumn(b)
        error("Errore. b = %s dovrebbe essere un vettore colonna.", stampaInLinea(b));
    end
    [nVincoli, ~] = size(A);
    if nVincoli ~= length(b)
        error("Errore. Il numero di vincoli (righe di A) dev'essere uguale alla lunghezza di b.");
    end
end

