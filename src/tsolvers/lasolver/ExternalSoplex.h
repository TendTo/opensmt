#ifndef EXTERNAL_SOPLEX_H
#define EXTERNAL_SOPLEX_H

#include "ExternalLpSolver.h"

#ifdef OPENSSMT_ENABLE_SOPLEX

namespace opensmt::externallp {

class ExternalSoplex : public ExternalLpSolver {
public:
    explicit ExternalSoplex(Simplex const & simplex, SMTConfig const & config);
    ~ExternalSoplex() override = default;

    Solution solve() override { return {}; }
};

class ExternalSoplexEpsilon final : public ExternalSoplex {
public:
    ExternalSoplexEpsilon(Simplex const & simplex, SMTConfig const & config);

protected:
    Solution extractSolution() override { return {}; }
};

class ExternalSoplexStrict final : public ExternalSoplex {
public:
    ExternalSoplexStrict(Simplex const & simplex, SMTConfig const & config);

protected:
    Solution extractSolution() override { return {}; }
};

} // namespace opensmt::externallp

#endif

namespace opensmt::externallp {
std::unique_ptr<ExternalLpSolver> make_external_soplex(Simplex const & simplex, SMTConfig const & config);
}

#endif
