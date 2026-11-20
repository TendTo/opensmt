#ifndef EXTERNAL_SOPLEX_H
#define EXTERNAL_SOPLEX_H

#define SOPLEX_WITH_GMP
#define SOPLEX_WITH_BOOST
#define SOPLEX_WITH_LIGHT_RAT

#include "ExternalLpSolver.h"
#include <soplex.h>

#ifdef OPENSSMT_ENABLE_SOPLEX

namespace opensmt::externallp {

class ExternalSoplex : public ExternalLpSolver {

public:
    using SpxStatus = soplex::SPxSolverBase<double>::Status;
    using VarStatus = soplex::SPxSolverBase<double>::VarStatus;

    explicit ExternalSoplex(Simplex const & simplex, SMTConfig const & config);

    Solution solve() override;

protected:
    [[nodiscard]] int numCols() const { return d_spx.numCols(); }
    [[nodiscard]] int numRows() const { return d_spx.numRows(); }

    LinResult d_sol_status;
    soplex::SoPlex d_spx;
    soplex::VectorRational d_x;
    soplex::VectorRational d_y;

public:
    enum class VariableType {
        ROW,
        COL,
    };

protected:
    template<VariableType VarType>
    void extractVarValue(int idx, Solution & sol, soplex::VectorRational const * values = nullptr);
    [[nodiscard]] virtual bool isStrictVarZero() const = 0;
};

class ExternalSoplexEpsilon final : public ExternalSoplex {
public:
    static constexpr bool isEpsilon = true;
    ExternalSoplexEpsilon(Simplex const & simplex, SMTConfig const & config);

protected:
    [[nodiscard]] bool isStrictVarZero() const override { return false; }

    [[nodiscard]] Solution extractSolution() override;

private:
    static constexpr double SMALL_FIXED_EPSILON = .000000001;
};

class ExternalSoplexStrict final : public ExternalSoplex {
public:
    static constexpr bool isEpsilon = false;
    ExternalSoplexStrict(Simplex const & simplex, SMTConfig const & config);

protected:
    [[nodiscard]] bool isStrictVarZero() const override {
        if (d_x.dim() == 0) return false;
        return d_spx.hasSol() ? d_x[numCols() - 1].is_zero() : false;
    }

    [[nodiscard]] Solution extractSolution() override;
};

} // namespace opensmt::externallp

#endif

namespace opensmt::externallp {
std::unique_ptr<ExternalLpSolver> make_external_soplex(Simplex const & simplex, SMTConfig const & config);
}

#endif
