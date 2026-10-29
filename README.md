# Exact Simplex and Integer Programming in MATLAB

A linear and integer programming solver written from scratch in MATLAB for the Operations Research course, working in exact rational arithmetic through the Symbolic Math Toolbox so that every tableau, basis and cut is computed without rounding. `symlinprog` mirrors the interface of MATLAB's `linprog` (`min c'x` with inequality, equality and bound constraints) and solves it with the primal simplex; `symintlinprog` solves the integer version with Gomory cuts applied through the dual simplex.

## Structure
- `PrimoLivello/`: input checks, constraint simplification, basis search, pretty printing.
- `SecondoLivello/`: primal and dual simplex, auxiliary problems for a feasible starting basis, bound handling, Gomory cuts.
- `TerzoLivello/`: the two user-facing functions `symlinprog(c, Aleq, bleq, Aeq, beq, l, u)` and `symintlinprog(...)`.

## How to run
MATLAB with the Symbolic Math Toolbox. Add the three folders to the path and call `[x, fval] = symlinprog(c, Aleq, bleq)`; see the header comments of each function (in Italian) for the full argument list.
