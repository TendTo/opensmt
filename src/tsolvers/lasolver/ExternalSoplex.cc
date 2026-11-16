#include "ExternalSoplex.h"

#ifdef OPENSSMT_ENABLE_SOPLEX

namespace opensmt::externallp {
std::unique_ptr<ExternalLpSolver> make_external_soplex(SMTConfig const & config) {
    if (config.epsilon_mode()) return std::make_unique<ExternalSoplexEpsilon>(config);
    return std::make_unique<ExternalSoplexStrict>(config);
}
ExternalSoplex::ExternalSoplex(SMTConfig const & config) : ExternalLpSolver(config) {}
} // namespace opensmt::externallp

#else

namespace opensmt::externallp {
std::unique_ptr<ExternalLpSolver> make_external_soplex(SMTConfig const &) {
    throw std::runtime_error("opensmt was not build with support for the soplex lp solver");
}
} // namespace opensmt::externallp

#endif
