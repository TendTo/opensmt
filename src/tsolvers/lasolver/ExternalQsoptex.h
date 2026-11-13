#ifndef EXTERNAL_QSOPTEX_H
#define EXTERNAL_QSOPTEX_H

#include "ExternalLpSolver.h"
#include "common/numbers/FastRational.h"

#ifdef OPENSSMT_ENABLE_QSOPTEX
extern "C" {
#include <qsopt_ex/QSopt_ex.h>
}

// These #defines from <qsopt_ex/QSopt_ex.h> cause problems for us
// because they mess with SoPlex's enums.
#undef OPTIMAL
#undef DUAL_INFEASIBLE

namespace opensmt::externallp {

class rat_class {
public:
    rat_class() : rat() { rat_init(rat); }
    explicit rat_class(FastRational const & val) : rat_class() { *this = val; }
    explicit rat_class(std::int32_t const val) : rat_class() { *this = val; }
    explicit rat_class(std::int64_t const val) : rat_class() { *this = val; }
    explicit rat_class(double const val) : rat_class() { *this = val; }
    ~rat_class() { rat_clear(rat); }

    rat_class(rat_class const & o);
    rat_class(rat_class && o) noexcept;
    rat_class & operator=(rat_class const & o);
    rat_class & operator=(rat_class && o) noexcept;

    rat_class & operator=(FastRational const & val);
    rat_class & operator=(std::int32_t val);
    rat_class & operator=(std::int64_t val);
    rat_class & operator=(double val);

    [[nodiscard]] int sgn() const { return rat_sgn(rat); }
    [[nodiscard]] rat_t const & data() const { return rat; }

    operator FastRational() const;

private:
    rat_t rat;
};

class ExternalQsoptex : public ExternalLpSolver {
public:
    explicit ExternalQsoptex(SMTConfig const & config);
    using VarStatus = char;
    using SolverStatus = int;
    ~ExternalQsoptex() override;

protected:
    int numCols() const { return rat_QSget_colcount(d_qsx); }
    int numRows() const { return rat_QSget_rowcount(d_qsx); }
#if 0

    mpq_class const & varToLb(ArithVar v) const;
    mpq_class const & varToUb(ArithVar v) const;
    bool hasStrictBound(ArithVar v) const;
    bool hasStrictLb(ArithVar v) const;
    bool hasStrictUB(ArithVar v) const;

    mpq_class const & getY(int const rowIdx) const {
        Assert(d_qsx && d_qsx->lp && d_qsx->lp->pIpiz);
        return toMpqClass(d_qsx->lp->pIpiz[rowIdx]);
    }
    virtual external::Solution extractSolution(bool mip) = 0;


    ArithVar _getArithVar(int nid, int M, int ind) const;
    ArithVar getArithVarFromRow(int nid, int ind) const {
        if (ind >= 0) {
            NodeLog const & nl = d_log.getNode(nid);
            return nl.lookupRowId(ind);
        }
        return ARITHVAR_SENTINEL;
    }

    ArithVar getArithVarFromStructural(int ind) const {
        if (ind >= 0) {
            unsigned u = (unsigned)ind;
            if (u < d_colToArithVar.size()) { return d_colToArithVar[u]; }
        }
        return ARITHVAR_SENTINEL;
    }
#endif

    rat_QSprob d_qsx;

    SolverStatus d_status;
    std::vector<rat_class> d_rhs;
    std::vector<rat_class> d_range;
    std::vector<char> d_sense;
    std::vector<rat_class> d_x;
    QSbasis d_basis;

public:
    enum class VariableType {
        ROW,
        COL,
    };

protected:
    template<VariableType VarType>
    Solution extractVarValue(int idx);

    virtual bool isStrictVarZero() const = 0;

public:
    void init() override {}
    Solution solve() override { throw std::runtime_error("ExternalQsoptex::solve() not implemented"); }

protected:
    bool d_solvedRelaxation;
    bool d_solvedMIP;
};

class ExternalQsoptexEpsilon : public ExternalQsoptex {
public:
    using ExternalQsoptex::ExternalQsoptex;

private:
    /** UTILITIES FOR DEALING WITH ESTIMATES */

    static constexpr double SMALL_FIXED_EPSILON = .000000001;

    bool isStrictVarZero() const override { return false; }
};

class ExternalQsoptexStrict : public ExternalQsoptex {
public:
    using ExternalQsoptex::ExternalQsoptex;

private:
    std::vector<LVRef> d_strictVars;

    bool isStrictVarZero() const override {
        if (!d_x.empty()) return d_x[numCols() - 1].sgn() == 0;
        return false;
    }
};
} // namespace opensmt::externallp

#endif

namespace opensmt::externallp {
std::unique_ptr<ExternalLpSolver> make_external_qsoptex(SMTConfig const & config);
}

#endif
