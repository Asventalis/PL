function [] = controlloNumSym(x)
    %CONTROLLONUMSYM Controlla se l'input è numerico o simbolico.
    %   Questa funzione controlla se l'input è una variabile numerica o una
    %   variabile simbolica di tipo numerico intero o razionale.
    %   PARAMETRI
    %   x: variabile da controllare
    
    if ~exist("x", "var")
        error("Errore. Serve un parametro alla funzione.");
    end
    if ~(isnumeric(x) || (class(x)=="sym" && all(isSymType(x,"real"),"all")))
        error("Errore. L'input deve essere un tipo numerico o un tipo simbolico rappresentante un numero razionale");
    end
end

