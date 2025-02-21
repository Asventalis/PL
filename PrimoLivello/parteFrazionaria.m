function [xFrazionaria] = parteFrazionaria(x)
% PARTEFRAZIONARIA Calcola la parte frazionaria di x.
%   La parte frazionaria è definita come x-floor(x) ed è sempre positiva,
%   dato che floor(x) è sempre più piccolo di x, essendo l'approssimazione
%   per difetto. Ciò vale sempre, anche per x negativi.
%   PARAMETRI
%   x: variabile di cui calcolare la parte frazionaria
%   OUTPUT
%   xFrazionaria: parte frazionaria di x
    if ~exist("x", "var")
        error("Errore. Serve un parametro alla funzione.");
    end
    controlloNumSym(x);
    xFrazionaria = x - floor(x);
end

