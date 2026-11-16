#include "ExternalLpSolver.h"
#include "ExternalQsoptex.h"
#include "ExternalSoplex.h"

#include <iostream>

namespace opensmt::externallp {

ExternalLpSolverStats ExternalLpSolver::stats;

std::unique_ptr<ExternalLpSolver> ExternalLpSolver::get_external_lp_solver(SMTConfig const & config) {
    switch (config.external_lp_solver()) {
        case opensmt::ExternalLpSolver::qsoptex:
            return make_external_qsoptex(config);
        case opensmt::ExternalLpSolver::soplex:
            return make_external_soplex(config);
        default:
            throw std::runtime_error("Unknown lp solver");
    }
}
ExternalLpSolver::ExternalLpSolver(SMTConfig const & config) : config(config), nnz(0) {}

void ExternalLpSolver::init(Simplex const & simplex) {
    Tableau const & tableau = simplex.getTableau();
    if (tableau.getNumOfCols() == 0 or tableau.getRows().empty()) return;

    col_to_var.reserve(tableau.numActiveColumns() - tableau.numActiveRows()); // No point in adding rows as columns
    row_to_var.reserve(tableau.numActiveRows());
    for (auto const & row : tableau.getRows()) {
        if (row) nnz += row->size();
    }

    var_to_col = std::vector<int>(tableau.getNumOfCols(), -1);
    for (std::uint32_t v = 0; v < tableau.getNumOfCols(); v++) {
        // Basic variables will form the rows of the problem, non-basic ones will be the
        if (tableau.isBasic(LVRef{.x = v})) {
            row_to_var.emplace_back(v);
        } else if (tableau.isNonBasic(LVRef{.x = v}) || tableau.isQuasiBasic(LVRef{.x = v})) {
            col_to_var.emplace_back(v);
        }
    }

    initImpl(simplex);

#if 0 && !defined(NDEBUG)
    tableau.print();
    std::cout << "Printing rows" << std::endl;
    for (auto const & row : tableau.getRows()) {
        if (row) row->print();
        if (row) row->size();
    }
    std::cout << "Printing cols" << std::endl;
    for (std::uint32_t i = 0; i < tableau.getNumOfCols(); ++i) {
        if (simplex.hasLBound(LVRef{.x = i})) {
            simplex.Lb(LVRef{.x = i}).print(std::cout);
        } else {
            std::cout << "-inf";
        }
        std::cout << " <= Var" << i << " <= ";
        if (simplex.hasUBound(LVRef{.x = i})) {
            simplex.Ub(LVRef{.x = i}).print(std::cout);
        } else {
            std::cout << "inf";
        }
        std::cout << std::endl;
    }
#endif
}
} // namespace opensmt::externallp
