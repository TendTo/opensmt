#include "ExternalLpSolver.h"
#include "ExternalQsoptex.h"
#include "ExternalSoplex.h"

#include <iostream>
#include <ranges>

namespace opensmt::externallp {

ExternalLpSolverStats ExternalLpSolver::stats;

std::unique_ptr<ExternalLpSolver> ExternalLpSolver::get_external_lp_solver(Simplex const & simplex,
                                                                           SMTConfig const & config) {
    switch (config.external_lp_solver()) {
        case opensmt::ExternalLpSolver::qsoptex:
            return make_external_qsoptex(simplex, config);
        case opensmt::ExternalLpSolver::soplex:
            return make_external_soplex(simplex, config);
        default:
            throw std::runtime_error("Unknown lp solver");
    }
}

void ExternalLpSolver::updateNbValue(Simplex & simplex, Solution const & sol, LVRef v) {
    if (sol.newBounds.contains(v)) {
        assert(sol.newBounds.at(v) == VarBoundState::VAR_BOUND_AT_LOWER ||
               sol.newBounds.at(v) == VarBoundState::VAR_BOUND_AT_UPPER);
        if (sol.newBounds.at(v) == VarBoundState::VAR_BOUND_AT_LOWER)
            simplex.changeValueToLb(v);
        else
            simplex.changeValueToUb(v);
        return;
    }
    assert(sol.newValues.contains(v));
    simplex.changeValueTo(v, sol.newValues.at(v));
}
bool ExternalLpSolver::isImpossibleToFix(Simplex const & simplex, LVRef const v) {
    assert(simplex.tableau.isBasic(v));
    auto const max_var_id = std::numeric_limits<unsigned>::max();
    LVRef y_found = LVRef::Undef;
    // Model doesn't fit the lower bound
    if (simplex.isModelOutOfLowerBound(v)) {
        // For the Bland rule
        auto curr_var_id_y = max_var_id;
        // look for nonbasic terms to fix the breaking of the bound
        for (auto const & term : simplex.tableau.getRowPoly(v)) {
            auto y = term.var;
            assert(v != y);
            assert(simplex.tableau.isNonBasic(y));
            auto const & coeff = term.coeff;
            bool const coeff_is_pos = isPositive(coeff);
            if ((coeff_is_pos && simplex.isModelStrictlyUnderUpperBound(y)) ||
                (!coeff_is_pos && simplex.isModelStrictlyOverLowerBound(y))) {
                // Choose the leftmost nonbasic variable with a negative (reduced) cost
                return false;
            }
        }
    } else if (simplex.isModelOutOfUpperBound(v)) {
        auto curr_var_id_y = max_var_id;
        // look for nonbasic terms to fix the unbounding
        for (auto const & term : simplex.tableau.getRowPoly(v)) {
            auto y = term.var;
            assert(v != y);
            assert(simplex.tableau.isNonBasic(y));
            auto const & coeff = term.coeff;
            bool const coeff_is_pos = isPositive(coeff);
            if ((!coeff_is_pos && simplex.isModelStrictlyUnderUpperBound(y)) ||
                (coeff_is_pos && simplex.isModelStrictlyOverLowerBound(y))) {
                return false;
            }
        }
    } else {
        assert(false);
    }
    return true;
}

/*
 * Make the currently basic variable v with a non basic variable in its row,
 * making v non-basic and the new variable basic.
 * - v must be a basic variable
 */
bool ExternalLpSolver::alignmentPivot(Simplex & simplex, Solution & sol, LVRef v) {
    auto const & tableau = simplex.getTableau();
    assert(tableau.isBasic(v));

    LVRef best_col = LVRef::Undef;
    std::uint32_t best_size = -1;

    // find the lower-cost variable to pivot
    auto const & row = tableau.getRowPoly(v);
    for (auto const & term : row) {
        assert(tableau.isNonBasic(term.var));
        // No point in making basic a variable that we want to be non-basic
        if (sol.newNonBasis.contains(term.var)) continue;
        auto const & col = tableau.getColumn(term.var);
        if (col.size() < best_size) {
            best_col = term.var;
            best_size = col.size();
        }
    }

    if (best_size == -1) return false;

    assert(!sol.newBounds.contains(v) || sol.newBounds.at(v) != VarBoundState::VAR_BOUND_AT_LOWER ||
           simplex.hasLBound(v));
    assert(!sol.newBounds.contains(v) || sol.newBounds.at(v) != VarBoundState::VAR_BOUND_AT_UPPER ||
           simplex.hasUBound(v));
    assert(sol.newBounds.contains(v) || (!simplex.hasLBound(v) && !simplex.hasUBound(v)) ||
           sol.newValues.at(v).R().isZero());

    // Determine the new value for the variable that is about to become non-basic
    Delta const & newVal =
        sol.newBounds.contains(v)
            ? (sol.newBounds.at(v) == VarBoundState::VAR_BOUND_AT_LOWER ? simplex.Lb(v) : simplex.Ub(v))
            : sol.newValues.at(v);
    simplex.pivot(v, best_col, &newVal);
    return true;
}

void ExternalLpSolver::importSolution(Simplex & simplex, Solution & sol) {
    // We only import solutions that are either feasible or infeasible, ignoring errors and similar situations
    if (sol.linResult != LinResult::LinFeasible && sol.linResult != LinResult::LinInfeasible) return;

    auto const & tableau = simplex.getTableau();

    for (std::uint32_t x = 0; x < tableau.getNumOfCols(); x++) {
        LVRef const v{.x = x};
        // Only update non-basic variables,
        // hoping to get a feasible solution immediately
        if (!tableau.isNonBasic(v)) continue;
        assert(tableau.isNonBasic(v));
        updateNbValue(simplex, sol, v);
        // The invariant could be broken temporarily due to importing a delta-complete model from a previous iteration,
        // but should hold after the importing
        assert(simplex.config.delta() >= 0 || simplex.invariantHolds());
    }

    if (simplex.config.delta() >= 0) {
        // If we are in delta-complete mode, some out-of-bounds basic variables may not be included in the candidates.
        // We need to sync that.
        for (std::uint32_t x = 0; x < tableau.getNumOfCols(); x++) {
            LVRef const v{.x = x};
            // Only update basic variables, ensuring that, if they are out of bounds, they are candidates
            if (!tableau.isBasic(v)) continue;
            assert(tableau.isBasic(v));
            if (simplex.isModelOutOfBounds(v)) simplex.candidates.emplace(v);
        }
    }
    assert(simplex.invariantHolds());
    assert(std::ranges::all_of(simplex.candidates, [&simplex](LVRef const var) {
        return var != LVRef::Undef && simplex.tableau.isBasic(var) && simplex.isModelOutOfBounds(var);
    }));
    assert(std::ranges::all_of(std::views::iota(0u, tableau.getNumOfCols()), [&simplex](std::uint32_t const var) {
        return !simplex.tableau.isBasic({.x = var}) || !simplex.isModelOutOfBounds({.x = var}) ||
               simplex.candidates.contains({.x = var});
    }));

    // No out-of-bounds basic variables, the problem is feasible
    if (simplex.candidates.empty()) { return; }
    // An out-of-bound basic variable cannot be fixed, the problem is infeasible.
    // Perform this expensive check only if we are confident the problem is infeasible
    if (sol.linResult == LinResult::LinInfeasible) {
        for (LVRef const v : simplex.candidates) {
            if (isImpossibleToFix(simplex, v)) {
                sol.infeasible_candidate = v;
                return;
            }
        }
    }

    // The simple assignment was not enough, extended search is needed
    stats.extended_search_count++;

    std::unordered_set<LVRef> missing_nb_vars;
    for (LVRef const v : sol.newNonBasis) {
        // Add the values that are basic but must become non-basic to the set
        if (tableau.isBasic(v)) missing_nb_vars.emplace(v);
    }

    bool performed_pivot = true;
    while (!missing_nb_vars.empty() && !simplex.candidates.empty() && performed_pivot) {
        performed_pivot = false;

        // Go over all variables that we need to add to the basis
        for (LVRef v : missing_nb_vars) {
            assert(tableau.isBasic(v));

            // Get the basic row
            if (alignmentPivot(simplex, sol, v)) {
                // A pivot has been performed.
                // v is now non-basic and can be removed from the set
                assert(tableau.isNonBasic(v));
                // An out-of-bound basic variable cannot be fixed, the problem is infeasible.
                // Perform this expensive check only if we are confident the problem is infeasible
                if (sol.linResult == LinResult::LinInfeasible) {
                    for (LVRef const c : simplex.candidates) {
                        if (isImpossibleToFix(simplex, c)) {
                            sol.infeasible_candidate = c;
                            return;
                        }
                    }
                }
                missing_nb_vars.erase(v);
                performed_pivot = true;
                break;
            }
        }

        assert(std::ranges::all_of(simplex.candidates, [&simplex](LVRef const var) {
            return var != LVRef::Undef && simplex.tableau.isBasic(var) && simplex.isModelOutOfBounds(var);
        }));
        assert(std::ranges::all_of(std::views::iota(0u, tableau.getNumOfCols()), [&simplex](std::uint32_t var) {
            return !simplex.tableau.isBasic({.x = var}) || !simplex.isModelOutOfBounds({.x = var}) ||
                   simplex.candidates.contains({.x = var});
        }));
        if (simplex.candidates.empty()) return;
    }

    assert(std::ranges::all_of(simplex.candidates, [&simplex](LVRef const var) {
        return var != LVRef::Undef && simplex.tableau.isBasic(var) && simplex.isModelOutOfBounds(var);
    }));
    assert(std::ranges::all_of(std::views::iota(0u, tableau.getNumOfCols()), [&simplex](std::uint32_t var) {
        return !simplex.tableau.isBasic({.x = var}) || !simplex.isModelOutOfBounds({.x = var}) ||
               simplex.candidates.contains({.x = var});
    }));
}

void ExternalLpSolver::copySolution(Simplex & simplex, Solution & sol) {
    // We only import solutions that are delta-feasible, ignoring errors and similar situations
    if (sol.linResult != LinResult::LinFeasible || simplex.config.delta() < 0) return;

    auto const & tableau = simplex.getTableau();

    for (std::uint32_t x = 0; x < tableau.getNumOfCols(); x++) {
        LVRef const v{.x = x};
        // Only update non-basic variables,
        // hoping to get a feasible solution immediately
        if (!tableau.isNonBasic(v)) continue;
        assert(tableau.isNonBasic(v));
        updateNbValue(simplex, sol, v);
    }
    // The invariant could be broken temporarily due to importing a delta-complete model from a previous iteration,
    // but should hold after the importing
    assert(simplex.invariantHolds());
}

ExternalLpSolver::ExternalLpSolver(Simplex const & simplex, SMTConfig const & config)
    : simplex(simplex),
      config(config),
      nnz(0) {
    Tableau const & tableau = simplex.getTableau();
    if (tableau.getNumOfCols() == 0 or tableau.getRows().empty()) return;

    col_to_var.reserve(tableau.numActiveColumns() - tableau.numActiveRows()); // No point in adding rows as columns
    row_to_var.reserve(tableau.numActiveRows());
    for (auto const & row : tableau.getRows()) {
        if (row) nnz += row->size();
    }

    var_to_col = std::vector<int>(tableau.getNumOfCols(), -1);
    for (std::uint32_t x = 0; x < tableau.getNumOfCols(); x++) {
        LVRef const v{.x = x};
        // Basic variables will form the rows of the problem, non-basic ones will be the
        if (tableau.isNonBasic(v)) {
            col_to_var.emplace_back(v);
        } else if (tableau.isBasic(v) && (simplex.hasLBound(v) || simplex.hasUBound(v))) { // Only add bounded rows
            row_to_var.emplace_back(v);
        }
        // Ignore quasi-basic variables, since they indicate deactivated constraints
    }

    stats.lp_max_nnz = std::max(stats.lp_max_nnz, nnz);
#if 1 && !defined(NDEBUG)
    tableau.print();
    // std::cout << "Printing rows" << std::endl;
    // for (auto const & row : tableau.getRows()) {
    //     if (row) row->print();
    // }
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
        std::cout << " -- "
                  << (tableau.isBasic({.x = i})      ? 'B'
                      : tableau.isNonBasic({.x = i}) ? 'N'
                                                     : 'Q')
                  << std::endl;
    }
#endif
}
#ifndef NDEBUG
std::ostream & operator<<(std::ostream & os, Solution const & sol) {
    os << "NewBounds: " << sol.newBounds << "\nCandidate: " << sol.infeasible_candidate.x
       << "\nNewNonBasis: " << sol.newNonBasis << "\nNewValues: " << sol.newValues << "\nLinResult: ";
    switch (sol.linResult) {
        case LinResult::LinError:
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
