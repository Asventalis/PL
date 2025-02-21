function [res] = simplessoPrimale(c,A,b,B,maxIter)
% SIMPLESSOPRIMALE Risolve un problema di programmazione lineare in forma
% primale standard.
%   Segue l'algoritmo del simplesso primale per risolvere il problema
%   di massimizzazione di c'*x con vincoli A*x<=b. Usa B come base di
%   partenza. B deve essere ammissibile per il primale. Se B non viene
%   passata, genera una base ammissibile per il primale con 
%   ausiliarioPrimale, se esiste. Ritorna res.val=NaN se non esiste
%   soluzione, res.val=inf se la soluzione è infinita, e ritorna un
%   valore finito se la soluzione esiste. Le altre info sono settate in
%   accordo. Ritorna inoltre steps, ovvero una struct contenente
%   informazioni per ogni passo dell'algoritmo che è stato svolto.
%   Steps può essere stampato in modo elegante sulla console
%   utilizzando stampaStruct.
%   PARAMETRI
%   c: vettore colonna della funzione obiettivo
%   A: matrice dei vincoli
%   b: vettore colonna dei termini noti
%   B: base ammissibile da cui partire
%   maxIter: numero massimo di iterazioni eseguibili
%   OUTPUT
%   res: struct:
%     - x: soluzione primale di base ottima
%     - y: soluzione duale di base ottima
%     - base: base ottima
%     - val: valore della soluzione di base ottima
%     - steps: struct:
%         - base: base dello step i-esimo
%         - invAB: l'inversa di A(base,:)
%         - x: soluzione primale della base
%         - y: soluzione duale della base
%         - h: indice uscente
%         - Wh: colonna di indice h della matrice W=-A(base,:)^-1
%         - AN_times_Wh: prodotto tra le righe non di base di A e Wh
%         - rapporti: i rapporti usati per calcolare k
%         - k: indice entrante
%   ESEMPIO
%   c = [4; 5; 2];
%   A = [0 0.6 0.8;
%       -1 2 0;
%       1 0 -1;
%       -1 0 0;
%       0 -1 0;
%       0 0 -1];
%   b = [500; 0; 0; 0; 0; 0];
%   B = [1 4 5];
%   res = simplessoPrimale(c,A,b,B);
%   stampaStruct(res.steps);
    if ~exist("c", "var") || ~exist("A", "var") || ~exist("b", "var") 
        error("Errore. c, A e b devono essere tutti e 3 specificati.");
    end
    controlloPrimale(c, A, b);
    [nVincoli,~] = size(A);
    if ~exist("B", "var") || isempty(B)
        B = ausiliarioPrimale(A, b);
        if isnan(B) % se ausiliarioPrimale torna NaN il problema è vuoto
            res.x = NaN;
            res.y = NaN;
            res.base = NaN;
            res.val = NaN;
            res.steps = NaN;
            return;
        end
    end
    controlloBase(B, A);
    if ~exist("maxIter","var")
        maxIter = -1;
    end
    if maxIter ~= -1 && (maxIter<0 || round(maxIter) ~= maxIter)
        error("Errore. MaxIter deve essere un numero naturale oppure -1.");
    end
    iter = 0;
    c = sym(c);
    A = sym(A);
    b = sym(b);
    x = A(B,:)^-1*b(B,:);
    if any(A*x-b>0)
        error("Errore. Base [%s] non ammissibile per il primale.", join(string(B),' '));
    end
    steps = [];
    step = struct;
    step.base = NaN;
    step.invAB = NaN;
    step.x = NaN;
    step.y = NaN;
    step.h = NaN;
    step.Wh = NaN;
    step.AN_times_Wh = NaN;
    step.rapporti = NaN;
    step.k = NaN;
    while maxIter == -1 || iter ~= maxIter
        iter = iter + 1;
        step.base = B;
        N = setdiff(1:nVincoli, B);
        step.invAB = A(B,:)^-1;
        W = -A(B,:)^-1;
        x = - W * b(B,:); % calcola x
        step.x = x;
        yB = -W' * c; % calcola yB
        y = sym(zeros(nVincoli, 1));
        y(B,:) = yB;
        step.y = y;
        if all(yB >= 0) % se yB ha tutti valori negativi, siamo all'ottimo
            steps = [steps; step];
            res.x = x;
            res.y = y;
            res.base = B;
            res.val = c' * x;
            res.steps = steps;
            return;
        end
        h = find(y<0,1); % altrimenti calcola h
        step.h = h;
        Wh = W(:,B==h); % calcola Wh
        step.Wh = Wh;
        step.AN_times_Wh = A(N,:)*Wh;
        if all(A(N,:)*Wh <= 0) % se tutti gli Ai*Wh con i in N sono negativi, la soluzione è infinito
            steps = [steps; step];
            Wh(Wh~=0) = Wh(Wh~=0)*inf;
            res.x = x + Wh;
            res.y = NaN;
            res.base = NaN;
            res.val = inf;
            res.steps = steps;
            return;
        end
        Npos = N(A(N,:)*Wh > 0);
        rapporti = (b(Npos,:)-A(Npos,:)*x)./(A(Npos,:)*Wh); % altrimenti calcola i rapporti
        step.rapporti = rapporti;
        rapportoMinimo = min(rapporti);
        k = Npos(find(rapporti==rapportoMinimo,1)); % calcola k
        step.k = k;
        steps = [steps; step];
        step = struct;
        step.base = NaN;
        step.invAB = NaN;
        step.x = NaN;
        step.y = NaN;
        step.h = NaN;
        step.Wh = NaN;
        step.AN_times_Wh = NaN;
        step.rapporti = NaN;
        step.k = NaN;
        B(B==h) = k; % metti k al posto di h
        B = sort(B);
    end
    res.x = x;
    res.y = y;
    res.base = B;
    res.val = c' * x;
    res.steps = steps;
end

