#include "ExternalSoplex.h"

#ifdef OPENSSMT_ENABLE_SOPLEX

namespace opensmt::externallp {
std::unique_ptr<ExternalLpSolver> make_external_soplex(Simplex const & simplex, SMTConfig const & config) {
    if (config.epsilon_mode()) return std::make_unique<ExternalSoplexEpsilon>(simplex, config);
    return std::make_unique<ExternalSoplexStrict>(simplex, config);
}
ExternalSoplex::ExternalSoplex(Simplex const & simplex, SMTConfig const & config) : ExternalLpSolver(simplex, config) {}
ExternalSoplexEpsilon::ExternalSoplexEpsilon(Simplex const & simplex, SMTConfig const & config)
    : ExternalSoplex(simplex, config) {}
ExternalSoplexStrict::ExternalSoplexStrict(Simplex const & simplex, SMTConfig const & config)
    : ExternalSoplex(simplex, config) {}
} // namespace opensmt::externallp

#else

namespace opensmt::externallp {
std::unique_ptr<ExternalLpSolver> make_external_soplex(SMTConfig const &) {
    throw std::runtime_error("opensmt was not build with support for the soplex lp solver");
}
} // namespace opensmt::externallp

#endif
