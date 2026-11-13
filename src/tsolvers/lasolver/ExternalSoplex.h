#ifndef EXTERNAL_SOPLEX_H
#define EXTERNAL_SOPLEX_H

#include "ExternalLpSolver.h"

#ifdef OPENSSMT_ENABLE_SOPLEX

namespace opensmt::externallp {

class ExternalSoplex : public ExternalLpSolver {
    using ExternalLpSolver::ExternalLpSolver;

    void init() override {}
    Solution solve() override { throw std::runtime_error("ExternalSoplex::solve() not implemented"); }
};

class ExternalSoplexEpsilon : public ExternalSoplex {
    using ExternalSoplex::ExternalSoplex;
};

class ExternalSoplexStrict : public ExternalSoplex {
    using ExternalSoplex::ExternalSoplex;
};

} // namespace opensmt::externallp

#endif

namespace opensmt::externallp {
std::unique_ptr<ExternalLpSolver> make_external_soplex(SMTConfig const & config);
}

#endif
