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
    void setInfinity();
    void setNInfinity();
    [[nodiscard]] bool isZero() const;
    [[nodiscard]] bool isNonZero() const;
    [[nodiscard]] bool isFinite() const;
    std::strong_ordering operator<=>(rat_class const & o) const;
    std::strong_ordering operator<=>(FastRational const & o) const;

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
    rat_t & operator*() { return rat; }
    rat_class operator+(rat_class const & o) const;
    rat_class operator-(rat_class const & o) const;
    rat_class operator-() const;

    operator FastRational() const;
    operator double() const;

private:
    rat_t rat;
};

#ifndef NDEBUG
std::ostream & operator<<(std::ostream & os, rat_class const & val);
std::ostream & operator<<(std::ostream & os, std::vector<rat_class> const & values);
#endif

class ExternalQsoptex : public ExternalLpSolver {
public:
    using VarStatus = char;
    using SolverStatus = int;

    explicit ExternalQsoptex(Simplex const & simplex, SMTConfig const & config);
    ~ExternalQsoptex() override;

    Solution solve() override;

protected:
    [[nodiscard]] int numCols() const { return rat_QSget_colcount(d_qsx); }
    [[nodiscard]] int numRows() const { return rat_QSget_rowcount(d_qsx); }
    virtual void solveImpl(unsigned int & precision, rat_class & delta) = 0;

    rat_QSprob d_qsx;

    SolverStatus d_status;
    std::vector<rat_class> d_rhs;
    std::vector<rat_class> d_range;
    std::vector<char> d_sense;
    std::vector<rat_class> d_x;
    std::vector<rat_class> d_y;
    QSbasis d_basis;

public:
    enum class VariableType {
        ROW,
        COL,
    };

protected:
    template<VariableType VarType>
    void extractVarValue(int idx, Solution & sol) const;
    [[nodiscard]] virtual bool isStrictVarZero() const = 0;
    [[nodiscard]] LinResult qsxStatusToSolStatus() const;
};

class ExternalQsoptexEpsilon final : public ExternalQsoptex {
public:
    static constexpr bool isEpsilon = true;
    ExternalQsoptexEpsilon(Simplex const & simplex, SMTConfig const & config);

protected:
    void solveImpl(unsigned int & precision, rat_class & delta) override;

    [[nodiscard]] bool isStrictVarZero() const override { return false; }

    [[nodiscard]] Solution extractSolution() override;

private:
    static constexpr double SMALL_FIXED_EPSILON = .000000001;
};

class ExternalQsoptexStrict final : public ExternalQsoptex {
public:
    static constexpr bool isEpsilon = false;
    ExternalQsoptexStrict(Simplex const & simplex, SMTConfig const & config);

protected:
    void solveImpl(unsigned int & precision, rat_class & delta) override;

    [[nodiscard]] bool isStrictVarZero() const override {
        if (!d_x.empty()) return d_x[numCols() - 1].sgn() == 0;
        return false;
    }

    [[nodiscard]] Solution extractSolution() override;
};
} // namespace opensmt::externallp

#endif

namespace opensmt::externallp {
std::unique_ptr<ExternalLpSolver> make_external_qsoptex(Simplex const & simplex, SMTConfig const & config);
}

#endif
