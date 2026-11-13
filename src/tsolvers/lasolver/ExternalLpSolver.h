#ifndef EXTERNAL_SIMPLEX_H
#define EXTERNAL_SIMPLEX_H

#include "LARefs.h"
#include "options/SMTConfig.h"

#include <memory>
#include <ostream>
#include <unordered_map>
#include <unordered_set>

namespace opensmt::externallp {

struct ExternalLpSolverStats {
    int external_pivots{0};
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
    int lp_max_nnz{0};
    int adjustment_pivots_feasible{0};
    int adjustment_pivots_infeasible{0};
    int lp_solver{0};
    int lp_setup_time{0};
    int lp_run_time{0};
    std::unordered_map<int, int> precisions{};
    std::unordered_map<int, int> iterative_refinements{};

    void printStatistics(std::ostream & os) const {
        os << "; -------------------------" << '\n';
        os << "; STATISTICS FOR EXT. SPX. " << '\n';
        os << "; -------------------------" << '\n';
        os << "; External pivots..........: " << external_pivots << '\n';
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
        os << "; LP setup time............: " << lp_setup_time << '\n';
        os << "; LP run time..............: " << lp_run_time << '\n';
        os << "; Precisions...............: {";
        for (auto const & [key, value] : precisions) {
            os << key << ": " << value << '\n';
        }
        os << "}\n";
        os << "; Iterative refinements....: {";
        for (auto const & [key, value] : iterative_refinements) {
            os << key << ": " << value << '\n';
        }
        os << "}\n";
        os << "; -------------------------" << '\n';
    }
};

enum LinResult {
    LinUnknown,    /* Unknown error */
    LinFeasible,   /* Relaxation is feasible */
    LinInfeasible, /* Relaxation is infeasible/all integer branches closed */
    LinExhausted
};

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
    std::vector<LVRef> newBasis;
    std::vector<LVRef> newNonBasis;
    std::vector<std::pair<LVRef, double>> newValues;
    LinResult linResult;
    Solution() : newBasis(), newNonBasis(), newValues(), linResult(LinUnknown) {}
};

class ExternalLpSolver {
public:
    static std::unique_ptr<ExternalLpSolver> get_external_lp_solver(SMTConfig const & config);

    explicit ExternalLpSolver(SMTConfig const & config) : config(config), stats() {}
    virtual ~ExternalLpSolver();

    virtual void init() = 0;
    virtual Solution solve() = 0;

protected:
    SMTConfig const & config;
    ExternalLpSolverStats stats;

    std::vector<LVRef> rowToArithVar;
    std::vector<LVRef> colToArithVar;
};

} // namespace opensmt::externallp

#endif
