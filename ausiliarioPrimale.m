function [baseAmm,info] = ausiliarioPrimale(A,b,B)
    %AUSILIARIOPRIMALE Trova una base ammissibile nel primale per un 
    %problema di programmazione lineare in forma primale standard.
    %   Trova una base ammissibile costruendo il problema primale
    %   ausiliario a partire dalla base B, e procedendo a cercare la
    %   soluzione. Se la soluzione ottimale ha valore diverso da 0, ritorna
    %   NaN al posto di baseAmm, altrimenti ricava la baseAmmissibile dalla
    %   soluzione ottimale del problema ausiliario. Ritorna inoltre info
    %   varie sullo svolgimento dell'algoritmo.
    %   PARAMETRI
    %   A: matrice dei vincoli
    %   b: vettore dei termini noti
    %   B: base da cui costruire il problema ausiliario
    %   OUTPUT
    %   baseAmm: una base ammissibile per il problema primale, oppure NaN
    %   info: struct:
    %     - cAux: funzione obiettivo del problema ausiliario
    %     - AAux: matrice dei vincoli del problema ausiliario
    %     - bAux: matrice dei termini noti del problema ausiliario
    %     - baseAux: base ammissibile da cui partire per risolvere il
    %       problema ausiliario con il simplesso primale
    %     - resSimplesso: output di simplessoPrimale applicato al problema
    %       primale ausiliario
    %   ESEMPIO
    %   c = [4; 5; 2];
    %   A = [0 0.6 0.8;
    %       -1 2 0;
    %       1 0 -1;
    %       -1 0 0;
    %       0 -1 0;
    %       0 0 -1];
    %   b = [500; 0; 0; 0; 0; 0];
    %   B = [1 2 6];
    %  [baseAmm,info] = ausiliarioPrimale(A,b,B);
    %   stampaStruct(info.resSimplesso.steps);
    %   baseAmm
    
    if ~exist("A", "var") || ~exist("b", "var")
        error("Errore. A e b devono essere tutti e 2 specificati.");
    end
    controlloAb(A,b);
    A = sym(A);
    b = sym(b);
    [nVincoli, nVariabili] = size(A);
    if ~exist("B", "var")
        B = trovaBase(A);
    end
    controlloBase(B,A);
    N = setdiff(1:nVincoli, B);
    x = A(B,:)^-1*b(B,:);
    if nVincoli == nVariabili || all(A*x-b<=0)
        info = NaN;
        baseAmm = B;
        return;
    end
    U = intersect(find((A*x-b)'<=0), N);
    V = intersect(find((A*x-b)'>0), N);
    nVarAux = length(V);
    AAux = [A zeros(nVincoli,nVarAux); zeros(nVarAux,nVariabili+nVarAux)];
    AAux(V,nVariabili+1:nVariabili+nVarAux) = -eye(nVarAux);
    AAux(nVincoli+1:nVincoli+nVarAux,nVariabili+1:nVariabili+nVarAux) = -eye(nVarAux);
    bAux = [b; zeros(nVarAux,1)];
    cAux = sym([zeros(nVariabili,1); -ones(nVarAux,1)]);
    baseAux = union(B, V);
    res = simplessoPrimale(cAux, AAux, bAux, baseAux); % se non iniziassimo da baseAux si creerebbero chiamate ricorsive infinite
    if res.val ~= 0 % se il valore ottimo è diverso da 0, non esiste base ammissibile per il primale
        baseAmm = NaN;
    else % se il valore ottimo è 0, si può trovare una base ammissibile a partire dalla base ottima del primale ausiliario
        baseAmpliata = res.base;
        % algoritmo per trovare una base ammissibile a partire da una base ottima dell'ausiliario (ipotizziamo che la dimostrazione nel libro del prof sia giusta)
        baseBU = intersect(baseAmpliata, union(B,U));
        baseV = intersect(baseAmpliata, V);
        baseAmm = union(baseBU, baseV(1:nVariabili-length(baseBU)));
    end
    info.cAux = cAux;
    info.AAux = AAux;
    info.bAux = bAux;
    info.baseAux = baseAux;
    info.resSimplesso = res;
end

