#include "ExternalLpSolver.h"
#include "ExternalQsoptex.h"
#include "ExternalSoplex.h"

#include <iostream>

namespace opensmt::externallp {
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
ExternalLpSolver::~ExternalLpSolver() {
#ifdef STATISTICS
    stats.printStatistics(std::cerr);
#endif // STATISTICS
}
} // namespace opensmt::externallp
