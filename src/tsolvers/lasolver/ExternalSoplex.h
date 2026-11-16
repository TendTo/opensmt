#ifndef EXTERNAL_SOPLEX_H
#define EXTERNAL_SOPLEX_H

#include "ExternalLpSolver.h"

#ifdef OPENSSMT_ENABLE_SOPLEX

namespace opensmt::externallp {

class ExternalSoplex : public ExternalLpSolver {
public:
    explicit ExternalSoplex(SMTConfig const & config);
    ~ExternalSoplex() override = default;

    Solution solve() override { return {}; }
};

class ExternalSoplexEpsilon final : public ExternalSoplex {
public:
    using ExternalSoplex::ExternalSoplex;

protected:
    void initImpl(Simplex const & simplex) override {}
};

class ExternalSoplexStrict final : public ExternalSoplex {
public:
    using ExternalSoplex::ExternalSoplex;

protected:
    void initImpl(Simplex const & simplex) override {}
};

} // namespace opensmt::externallp

#endif

namespace opensmt::externallp {
std::unique_ptr<ExternalLpSolver> make_external_soplex(SMTConfig const & config);
}

#endif
