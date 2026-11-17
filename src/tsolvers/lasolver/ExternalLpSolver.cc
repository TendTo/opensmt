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
#ifndef NDEBUG
std::ostream & operator<<(std::ostream & os, Solution const & sol) {
    os << "NewBounds: " << sol.newBounds << "\nNewBasis: " << sol.newBasis << "\nNewNonBasis: " << sol.newNonBasis
       << "\nNewValues: " << sol.newValues << "\nLinResult: ";
    switch (sol.linResult) {
        case LinResult::LinUnknown:
            return os << "unknown";
        case LinResult::LinFeasible:
            return os << "feasible";
        case LinResult::LinInfeasible:
            return os << "infeasible";
        case LinResult::LinExhausted:
            return os << "exhausted";
        default:
            return os;
    }
}
std::ostream & operator<<(std::ostream & os, std::unordered_set<LVRef> const & set) {
    os << '{';
    bool first = true;
    for (LVRef const v : set) {
        if (!first) os << ", ";
        os << 'v' << v.x;
        first = false;
    }
    return os << '}';
}
std::ostream & operator<<(std::ostream & os, std::unordered_map<LVRef, Delta> const & map) {
    os << '{';
    bool first = true;
    for (auto const & [v, value] : map) {
        if (!first) os << ", ";
        os << 'v' << v.x << ": " << value;
        first = false;
    }
    return os << '}';
}
std::ostream & operator<<(std::ostream & os, std::unordered_map<LVRef, VarBoundState> const & map) {
    os << '{';
    bool first = true;
    for (auto const & [v, bound] : map) {
        if (!first) os << ", ";
        os << 'v' << v.x << ": ";
        switch (bound) {
            case VarBoundState::VAR_BOUND_AT_LOWER:
                os << "lower";
                break;
            case VarBoundState::VAR_BOUND_AT_UPPER:
                os << "upper";
                break;
        }
        first = false;
    }
    return os << '}';
}
std::ostream & operator<<(std::ostream & os, std::vector<int> values) {
    os << '[';
    bool first = true;
    for (int const value : values) {
        if (!first) os << ", ";
        os << value;
        first = false;
    }
    return os << ']';
}
std::ostream & operator<<(std::ostream & os, std::vector<char> values) {
    os << '[';
    bool first = true;
    for (char const value : values) {
        if (!first) os << ", ";
        os << value;
        first = false;
    }
    return os << ']';
}
#endif

} // namespace opensmt::externallp
