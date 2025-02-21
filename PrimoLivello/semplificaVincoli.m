function [ASemplificati,bSemplificati] = semplificaVincoli(A,b)
% SEMPLIFICAVINCOLI Semplifica i vincoli in modo da rendere tutti i
% coefficienti numeri interi.
%   Dati dei vincoli nella forma A*x<=b,A*x=b o A*x>=b, semplifica le
%   disequazioni o equazioni in modo che tutti i coefficienti, ovvero
%   gli elementi di A e b, siano numeri interi, senza cambiare gli
%   effettivi vincoli.
%   PARAMETRI
%   A: matrice dei vincoli
%   b: vettore dei termini noti
%   OUTPUT
%   ASemplificati: matrice dei vincoli semplificati
%   bSemplificati: vettore dei termini noti semplificati
    if ~exist("A","var") || ~exist("b","var")
        error("Errore. A e b devono essere tutti e 2 specificati.");
    end
    controlloNumSym(A);
    controlloNumSym(b);
    nVincoli = size(A,1);

    if nVincoli ~= length(b)
        error("Errore. Il numero di righe di A deve essere uguale alla lunghezza di b.");
    end
    if ~iscolumn(b)
        error("Errore. b deve essere un vettore colonna.");
    end

    A = sym(A);
    b = sym(b);
    commonFactor = zeros(nVincoli,1);
    for i = 1:nVincoli
    	commonFactor(i) = gcd([A(i,:) b(i)])^-1; % gcd ritorna sempre valori positivi, quindi non cambia il segno dell'uguaglianza. Vedere: https://math.stackexchange.com/questions/151081/gcd-of-rationals
    end
    ASemplificati = A .* commonFactor;
    bSemplificati = b .* commonFactor;
end

