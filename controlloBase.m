function [] = controlloBase(B,A)
    %CONTROLLOBASE Controlla se B è una possibile base per un problema di
    %programmazione lineare in forma primale standard con matrice dei
    %vincoli A.
    %   La matrice A può essere omessa per limitarsi a controllare che B
    %   sia una generica possibile base. Controlla le seguenti cose:
    %    - B è un vettore di numeri o simboli numerici
    %    - B è un vettore riga
    %    - B è un vettore di numeri naturali
    %    - B non contiene doppioni
    %   Se la matrice A è passata come argomento, allora vengono
    %   controllate anche le seguenti cose:
    %    - gli elementi di B sono compresi tra 1 e il numero di vincoli
    %    - il numero di elementi di B è esattamente il numero di variabili
    %    - la matrice di base corrispondente è invertibile
    %   PARAMETRI
    %   B: base di un problema di programmazione lineare
    %   A: matrice dei vincoli
    
    if ~exist("B", "var")
        error("Errore. La base B va passata come primo argomento.");
    end
    if isempty(B)
        error("Errore. La base B è vuota.");
    end
    controlloNumSym(B);
    if ~isrow(B)
        error("Errore. B deve essere un vettore riga.");
    end
    if ~all(round(B) == B & B>0)
        error("Errore. Gli elementi di B devono essere numeri naturali.");
    end
    if ~isequal(unique(B), B)
        error("Errore. Non possono esserci doppioni in B e gli indici devono essere in ordine crescente.");
    end
    if ~exist("A", "var")
        return;
    end
    controlloA(A);
    [nVincoli,nVariabili] = size(A);
    if ~all(1<=B | B<=nVincoli)
        error("Errore. I valori di B dovrebbero essere tra 1 e il numero di vincoli.");
    end
    if length(B) ~= nVariabili
        error("Errore. B dovrebbe avere un numero di elementi pari al numero di variabili.");
    end
    if rank(A(B,:)) ~= nVariabili
        error("Errore. La matrice di base A(B,:) non è invertibile.");
    end
end

