#include "ExternalQsoptex.h"

#include <utility>

#ifdef OPENSSMT_ENABLE_QSOPTEX

namespace opensmt::externallp {

namespace {

    constexpr int pivot_threshold = 1000000;

    rat_class minus_one;
    rat_class infinity;
    rat_class ninfinity;

    void QSXStart() {
        static bool started = false;
        if (started) return;
        started = true;
        QSexactStart();
        minus_one = -1;
        infinity.setInfinity();
        ninfinity.setNInfinity();
    }

    void QSXFinish() {
        QSexactClear();
    }

    template<bool epsilon_mode, bool is_lower_bound>
    rat_class const & get_bound(Simplex const & simplex, LVRef const v, double const delta, rat_class & bound) {
        if (bool const is_finite = is_lower_bound ? simplex.hasLBound(v) : simplex.hasUBound(v); !is_finite)
            return is_lower_bound ? infinity : ninfinity;

        assert(is_lower_bound ? simplex.hasLBound(v) : simplex.hasUBound(v));
        Delta const & delta_bound = is_lower_bound ? simplex.Lb(v) : simplex.Ub(v);
        bool const is_strict = delta < 0 && delta_bound.hasDelta();

        bound = delta_bound.R();
        if (epsilon_mode && is_strict) {
            // Adjust the bounds with epsilon
            if (is_lower_bound) {
                rat_add_si(*bound, ExternalLpSolver::epsilon_num, ExternalLpSolver::epsilon_den);
            } else {
                rat_sub_si(*bound, ExternalLpSolver::epsilon_num, ExternalLpSolver::epsilon_den);
            }
        }
        return bound;
    }

} // namespace

rat_class::rat_class(rat_class const & o) : rat_class() {
    *this = o;
}
rat_class::rat_class(rat_class && o) noexcept : rat_class() {
    *this = std::move(o);
}
rat_class & rat_class::operator=(rat_class const & o) {
    rat_set(rat, o.data());
    return *this;
}
rat_class & rat_class::operator=(rat_class && o) noexcept {
    std::swap(rat, o.rat);
    return *this;
}
rat_class & rat_class::operator=(FastRational const & val) {
    if (auto const pair = val.tryGetNumDen(); pair.has_value()) {
        rat_set_int32(rat, pair.value().first, pair.value().second);
    } else {
        rat_set_mpq(rat, val.getMpq().get_mpq_t());
    }
    return *this;
}
rat_class & rat_class::operator=(std::int32_t const val) {
    rat_set_int32(rat, val, 1);
    return *this;
}
rat_class & rat_class::operator=(std::int64_t const val) {
    rat_set_int64(rat, val, 1);
    return *this;
}
rat_class & rat_class::operator=(double const val) {
    rat_set_double(rat, val);
    return *this;
}
rat_class rat_class::operator+(rat_class const & o) const {
    rat_class res;
    rat_add(res.rat, rat, o.rat);
    return res;
}
rat_class rat_class::operator-(rat_class const & o) const {
    rat_class res;
    rat_sub(res.rat, rat, o.rat);
    return res;
}
rat_class rat_class::operator-() const {
    rat_class res = *this;
    rat_neg_self(res.rat);
    return res;
}
rat_class::operator FastRational() const {
    return rat_is_gmp(rat) ? FastRational(static_cast<mpq_class>(rat_get_mpq(rat)))
                           : FastRational{rat_get_num(rat), rat_get_den(rat)};
}
void rat_class::setInfinity() {
    rat_set_infinity(rat);
}
void rat_class::setNInfinity() {
    rat_set_ninfinity(rat);
}
bool rat_class::isZero() const {
    return rat_is_zero(rat);
}
bool rat_class::isNonZero() const {
    return !rat_is_zero(rat);
}

#ifndef NDEBUG
std::ostream & operator<<(std::ostream & os, rat_class const & val) {
    char * const str = rat_get_str(val.data());
    os << str;
    free(str);
    return os;
}
std::ostream & operator<<(std::ostream & os, std::vector<rat_class> const & values) {
    os << '[';
    bool first = true;
    for (rat_class const & value : values) {
        if (!first) os << ", ";
        os << value;
        first = false;
    }
    return os << ']';
}
#endif

ExternalQsoptex::ExternalQsoptex(Simplex const & simplex, SMTConfig const & config)
    : ExternalLpSolver(simplex, config),
      d_status(),
      d_basis{0} {
    QSXStart();
    d_qsx = rat_QScreate_prob(nullptr, QS_MIN);

    rat_QSset_param(d_qsx, QS_PARAM_SIMPLEX_MAX_ITERATIONS, pivot_threshold);
    rat_QSset_param(d_qsx, QS_PARAM_SIMPLEX_DISPLAY, config.verbosity());
}

ExternalQsoptex::~ExternalQsoptex() {
    rat_QSfree_prob(d_qsx);
    if (d_basis.cstat != nullptr) free(d_basis.cstat);
    if (d_basis.rstat != nullptr) free(d_basis.rstat);
    d_basis.cstat = nullptr;
    d_basis.rstat = nullptr;
    d_basis.nstruct = 0;
    d_basis.nrows = 0;
}
Solution ExternalQsoptex::solve() {
#if 1 && !defined(NDEBUG)
    static int i = 0;
    // rat_QSdump_prob(qsx.prob);
    fprintf(stderr, "Problem %d %s: n_cols = %d, n_rows = %d\n", i, config.epsilon_mode() ? "epsilon" : "strict",
            numCols(), numRows());
    char filename[256];
    snprintf(filename, sizeof(filename),
             "/home/campus.ncl.ac.uk/c3054737/Programming/phd/opensmt/"
             "qsproblem_%d_%s.lp",
             i, config.epsilon_mode() ? "epsilon" : "strict");
    rat_QSwrite_prob(d_qsx, filename, "LP");
    i++;
#endif
    return {};
}

void ExternalQsoptexEpsilon::initImpl(Simplex const & simplex) {
    Tableau const & tableau = simplex.getTableau();

    std::vector<int> numNonZeroPerRow;
    std::vector<int> beginRowIdx;
    std::vector<int> colIdxs;
    std::vector<rat_class> values;

    d_rhs.reserve(row_to_var.size());
    d_range.reserve(row_to_var.size());
    numNonZeroPerRow.reserve(row_to_var.size());
    beginRowIdx.reserve(row_to_var.size());
    colIdxs.reserve(nnz);
    values.reserve(nnz);
    d_sense.reserve(row_to_var.size());

    rat_class lb_temp, ub_temp;
    for (LVRef const v : col_to_var) {
        assert(tableau.isNonBasic(v) || tableau.isQuasiBasic(v));
        rat_class const & lb = get_bound<true, true>(simplex, v, config.delta(), lb_temp);
        rat_class const & ub = get_bound<true, false>(simplex, v, config.delta(), ub_temp);

        var_to_col[v.x] = numCols();
        rat_QSnew_col(d_qsx, rat_zeroLpNum, lb.data(), ub.data(), nullptr);
    }

    // Construct the rows of the LP by parsing the polynomial constraints together
    // with the row bounds on the non-basic variables
    for (LVRef const v : row_to_var) {
        assert(tableau.isBasic(v));

        PolynomialT<LVRef> const & p = tableau.getRowPoly(v);
        numNonZeroPerRow.emplace_back(p.size());
        beginRowIdx.emplace_back(colIdxs.size());
        for (auto const & term : p) {
            int const colIdx = var_to_col[term.var.x];
            assert(colIdx >= 0);
            colIdxs.emplace_back(colIdx);
            values.emplace_back(term.coeff);
        }

        rat_class const & lb = get_bound<true, true>(simplex, v, config.delta(), lb_temp);
        rat_class const & ub = get_bound<true, false>(simplex, v, config.delta(), ub_temp);
        if (simplex.hasLBound(v) && simplex.hasUBound(v)) // Bounded
        {
            d_rhs.emplace_back(lb);
            d_range.emplace_back(ub - lb);
            d_sense.emplace_back('R');
        } else if (simplex.hasLBound(v)) // Only lower bound
        {
            d_rhs.emplace_back(lb);
            d_range.emplace_back();
            d_sense.emplace_back('G');
        } else if (simplex.hasUBound(v)) // Only upper bound
        {
            d_rhs.emplace_back(ub);
            d_range.emplace_back();
            d_sense.emplace_back('L');
        } else // Unbounded
        {
            d_rhs.emplace_back(ninfinity);
            d_range.emplace_back();
            d_sense.emplace_back('G');
        }
    }

    static_assert(sizeof(rat_class) == sizeof(rat_t), "mpq_class layout assumption broken");
    rat_QSadd_ranged_rows(d_qsx, static_cast<int>(numNonZeroPerRow.size()), numNonZeroPerRow.data(), beginRowIdx.data(),
                          colIdxs.data(), reinterpret_cast<rat_t const *>(values.data()),
                          reinterpret_cast<rat_t const *>(d_rhs.data()), d_sense.data(),
                          reinterpret_cast<rat_t const *>(d_range.data()), nullptr);
}

std::unique_ptr<ExternalLpSolver> make_external_qsoptex(SMTConfig const & config) {
    if (config.epsilon_mode()) return std::make_unique<ExternalQsoptexEpsilon>(config);
    return std::make_unique<ExternalQsoptexStrict>(config);
}

} // namespace opensmt::externallp

#else

namespace opensmt::externallp {
std::unique_ptr<ExternalLpSolver> make_external_qsoptex(SMTConfig const &) {
    throw std::runtime_error("opensmt was not build with support for the qsoptex lp solver");
}
} // namespace opensmt::externallp

#endif
