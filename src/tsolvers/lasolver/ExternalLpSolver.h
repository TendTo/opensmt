#ifndef EXTERNAL_SIMPLEX_H
#define EXTERNAL_SIMPLEX_H

#include "LARefs.h"
#include "Simplex.h"
#include "options/SMTConfig.h"

#include <memory>
#include <ostream>
#include <unordered_map>
#include <unordered_set>

namespace opensmt::externallp {

struct ExternalLpSolverStats {
    double max_delta{0.0};
    int delta_calls{0};
    int calls{0};
    int extended_search_count{0};
    int exhausted{0};
    int feasible_correct{0};
    int feasible_incorrect{0};
    int infeasible_correct{0};
    int infeasible_incorrect{0};
    int other{0};
    std::size_t lp_max_nnz{0};
    int adjustment_pivots_feasible{0};
    int adjustment_pivots_infeasible{0};
    int lp_solver{0};
    double lp_setup_time{0};
    double lp_run_time{0};
    std::unordered_map<unsigned int, int> precisions{};
    std::unordered_map<int, int> iterative_refinements{};

    void printStatistics(std::ostream & os) const {
        if (calls <= 0) return;
        os << "; -------------------------" << '\n';
        os << "; STATISTICS FOR EXT. SPX. " << '\n';
        os << "; -------------------------" << '\n';
        os << "; Max delta................: " << max_delta << '\n';
        os << "; Delta calls..............: " << delta_calls << '\n';
        os << "; External calls...........: " << calls << '\n';
        os << "; Extended search count....: " << extended_search_count << '\n';
        os << "; Exhausted................: " << exhausted << '\n';
        os << "; Feasible correct.........: " << feasible_correct << '\n';
        os << "; Feasible incorrect.......: " << feasible_incorrect << '\n';
        os << "; Infeasible correct.......: " << infeasible_correct << '\n';
        os << "; Infeasible incorrect.....: " << infeasible_incorrect << '\n';
        os << "; Other....................: " << other << '\n';
        os << "; LP max nnz...............: " << lp_max_nnz << '\n';
        os << "; LP solver................: " << lp_solver << '\n';
        os << "; LP setup time............: " << lp_setup_time << " s\n";
        os << "; LP run time..............: " << lp_run_time << " s \n";
        os << "; Adjustment pivots feas...: " << adjustment_pivots_feasible << '\n';
        os << "; Adjustment pivots infeas.: " << adjustment_pivots_infeasible << '\n';
        os << "; Precisions...............: {";
        for (auto const & [key, value] : precisions) {
            os << key << ": " << value << ", ";
        }
        os << "}\n";
        os << "; Iterative refinements....: {";
        for (auto const & [key, value] : iterative_refinements) {
            os << key << ": " << value << ", ";
        }
        os << "}\n";
        os << "; -------------------------" << '\n';
    }

    ~ExternalLpSolverStats() { printStatistics(std::cerr); }
};

enum class LinResult {
    LinUnknown,    /* The solver has not been invoked yet */
    LinEmtpy,      /* There are no constraints */
    LinError,      /* An error has occurred */
    LinFeasible,   /* Relaxation is feasible */
    LinInfeasible, /* Relaxation is infeasible/all integer branches closed */
    LinExhausted
};

enum class VarBoundState { VAR_BOUND_AT_LOWER, VAR_BOUND_AT_UPPER };

/**
 * Store a solution obtained from the external simplex solver.
 * It will contain all the information needed to attempt to apply the solution
 * to the current tableau, and to extract a model if the solution is integer
 * feasible.
 * Note that only one between `newBasis` and `newNonBasis` will be populated,
 * depending on the heuristic the external simplex solver wants the internal
 * tableau to use.
 * Moreover, only one between `linResult` and `mipResult` may differ from
 * `Unknown`, depending on the type of solution obtained from the external
 * simplex solver.
 */
struct Solution {
    // std::unordered_set<LVRef> newBasis;
    std::unordered_set<LVRef> newNonBasis;
    std::unordered_map<LVRef, Delta> newValues;
    std::unordered_map<LVRef, VarBoundState> newBounds;
    LinResult linResult;
    LVRef infeasible_candidate;
    Solution(LinResult const result = LinResult::LinUnknown, int const num_cols = 0, int const num_rows = 0)
        : linResult(result),
          infeasible_candidate(LVRef::Undef) {
        newNonBasis.reserve(num_cols);
        newValues.reserve(num_cols + num_rows);
        newBounds.reserve(num_cols + num_rows);
    }
};

class ExternalLpSolver {
public:
    static constexpr int lp_pivot_threshold = 1000000;
    static constexpr std::int32_t epsilon_num = 1;
    static constexpr std::uint32_t epsilon_den = 1000000000;
    static ExternalLpSolverStats stats;
    static std::unique_ptr<ExternalLpSolver> get_external_lp_solver(Simplex const & simplex, SMTConfig const & config);
    static void updateNbValue(Simplex & simplex, Solution const & sol, LVRef v);
    static bool isImpossibleToFix(Simplex const & simplex, LVRef v);
    static void importSolution(Simplex & simplex, Solution & sol);
    static bool alignmentPivot(Simplex & simplex, Solution & sol, LVRef v);

    explicit ExternalLpSolver(Simplex const & simplex, SMTConfig const & config);
    virtual ~ExternalLpSolver() = default;

    [[nodiscard]] bool hasBothBounds(LVRef const v) const { return simplex.hasLBound(v) && simplex.hasUBound(v); }
    [[nodiscard]] bool hasEitherBound(LVRef const v) const { return simplex.hasLBound(v) || simplex.hasUBound(v); }
    [[nodiscard]] bool hasStrictLb(LVRef const v) const { return simplex.hasLBound(v) && !simplex.Lb(v).D().isZero(); }
    [[nodiscard]] bool hasStrictUb(LVRef const v) const { return simplex.hasUBound(v) && !simplex.Ub(v).D().isZero(); }
    [[nodiscard]] bool hasStrictBound(LVRef const v) const { return hasStrictLb(v) || hasStrictUb(v); }
    [[nodiscard]] bool hasBothStrictBound(LVRef const v) const { return hasStrictLb(v) && hasStrictUb(v); }

    virtual Solution solve() = 0;

protected:
    virtual Solution extractSolution() = 0;

    Simplex const & simplex;
    SMTConfig const & config;

    std::vector<LVRef> row_to_var;
    std::vector<LVRef> col_to_var;
    std::vector<int> var_to_col;
    // std::vector<int> var_to_row;
    std::size_t nnz;
};

#ifndef NDEBUG

std::ostream & operator<<(std::ostream & os, Solution const & sol);
std::ostream & operator<<(std::ostream & os, std::unordered_set<LVRef> const & set);
std::ostream & operator<<(std::ostream & os, std::unordered_map<LVRef, Delta> const & map);
std::ostream & operator<<(std::ostream & os, std::unordered_map<LVRef, VarBoundState> const & map);
std::ostream & operator<<(std::ostream & os, std::vector<int>);
std::ostream & operator<<(std::ostream & os, std::vector<char>);

#endif

} // namespace opensmt::externallp

#endif
