#include "ExternalQsoptex.h"

#include <utility>

#ifdef OPENSSMT_ENABLE_QSOPTEX

namespace opensmt::externallp {

namespace {

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
            return is_lower_bound ? ninfinity : infinity;

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

    void ensure_perturbed_bounds_do_not_overlap(ExternalLpSolver const & solver, LVRef v, rat_class const & lb,
                                                rat_class const & ub, rat_class & lb_temp, rat_class & ub_temp) {
        if (solver.hasStrictBound(v) && lb > ub) {
            assert(lb.isFinite() && ub.isFinite());
            if (!solver.hasStrictLb(v)) {
                assert(&ub == &ub_temp);
                ub_temp = lb;
            } else {
                assert(&lb == &lb_temp);
                lb_temp = ub;
            }
        }
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
rat_class::operator double() const {
    return rat_to_d(rat);
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
bool rat_class::isFinite() const {
    return rat_is_finite(rat);
}

std::strong_ordering rat_class::operator<=>(rat_class const & o) const {
    int const res = rat_cmp(rat, o.rat);
    if (res == 0) return std::strong_ordering::equal;
    return res < 0 ? std::strong_ordering::less : std::strong_ordering::greater;
}

std::strong_ordering rat_class::operator<=>(FastRational const & o) const {
    int res;
    if (auto const pair = o.tryGetNumDen(); pair.has_value()) {
        res = rat_cmp_int32(rat, pair.value().first, pair.value().second);
    } else {
        res = rat_cmp_mpq(rat, o.getMpq().get_mpq_t());
    }
    if (res == 0) return std::strong_ordering::equal;
    return res < 0 ? std::strong_ordering::less : std::strong_ordering::greater;
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
std::ostream & operator<<(std::ostream & os, std::vector<mpq_class> const & values) {
    os << '[';
    bool first = true;
    for (mpq_class const & value : values) {
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

    rat_QSset_param(d_qsx, QS_PARAM_SIMPLEX_MAX_ITERATIONS, lp_pivot_threshold);
    rat_QSset_param(d_qsx, QS_PARAM_SIMPLEX_DISPLAY, config.verbosity());

    d_y.resize(row_to_var.size());
    d_x.resize(col_to_var.size() + 1); // Leave space for the strict-var
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
    if (col_to_var.empty() || row_to_var.empty()) return {LinResult::LinEmtpy};

    unsigned int precision = 0;
    rat_class delta{std::max(config.delta(), 0.0)};
    solveImpl(precision, delta);

    // Update the statistics
    stats.precisions[precision]++;
    if (config.delta() >= 0) {
        stats.max_delta = std::max(stats.max_delta, static_cast<double>(delta));
        if (delta.isNonZero()) stats.delta_calls++;
    }
    return extractSolution();
}

template<ExternalQsoptex::VariableType VarType>
void ExternalQsoptex::extractVarValue(int const idx, Solution & sol) const {
    auto & newNonBasis = sol.newNonBasis;
    auto & newValues = sol.newValues;
    auto & newBounds = sol.newBounds;

    LVRef v = LVRef::Undef;
    VarStatus varStatus = 'x';

    if constexpr (VarType == VariableType::COL) {
        v = col_to_var.at(idx);
        varStatus = d_basis.cstat[idx];
    }
    if constexpr (VarType == VariableType::ROW) {
        v = row_to_var.at(idx);
        varStatus = d_basis.rstat[idx];
    }
    assert(v != LVRef::Undef);
    assert(varStatus != 'x');

    if (varStatus != QS_COL_BSTAT_BASIC) newNonBasis.insert(v);

    rat_class const * value = nullptr;
    rat_class rhs;
    switch (varStatus) {
        // If we are dealing with a basic variable, we necessarily need to get its
        // value from the solved problem.
        case QS_COL_BSTAT_BASIC:
        case QS_COL_BSTAT_FREE: // Shared with basic
            if constexpr (VarType == VariableType::ROW) {
                rhs = d_rhs[idx] + (d_sense[idx] == 'L' ? -d_y[idx] : d_y[idx]);
                value = &rhs;
            } else if constexpr (VarType == VariableType::COL) {
                value = &d_x[idx];
            }
            assert(value != nullptr);
            if (simplex.hasLBound(v) && *value <= simplex.Lb(v).R()) {
                newBounds.emplace(v, VarBoundState::VAR_BOUND_AT_LOWER);
            } else if (simplex.hasUBound(v) && *value >= simplex.Ub(v).R()) {
                newBounds.emplace(v, VarBoundState::VAR_BOUND_AT_UPPER);
            } else {
                newValues.emplace(v, Delta(*value));
            }
            break;
        // For non-basic variables we know the value is at a bound
        // and we can use the d_vars/proper LP bounds directly
        case QS_COL_BSTAT_LOWER:
            // Free rows are set to lower bound, handle them accordingly
            if (!hasEitherBound(v)) {
                if (VarType == VariableType::ROW) {
                    rhs = d_rhs[idx] + (d_sense[idx] == 'G' ? d_y[idx] : -d_y[idx]);
                    value = &rhs;
                } else if (VarType == VariableType::COL) {
                    value = &d_x[idx];
                }
                assert(value != nullptr);
                newValues.emplace(v, Delta(*value));
            } else if constexpr (VarType == VariableType::ROW) {
                if (d_sense[idx] == 'G' || d_sense[idx] == 'R') {
                    assert(simplex.hasLBound(v));
                    newBounds.emplace(v, VarBoundState::VAR_BOUND_AT_LOWER);
                } else if (d_sense[idx] == 'L') {
                    assert(simplex.hasUBound(v));
                    newBounds.emplace(v, VarBoundState::VAR_BOUND_AT_UPPER);
                } else {
                    assert(false);
                }
            } else if constexpr (VarType == VariableType::COL) {
                assert(simplex.hasLBound(v));
                newBounds.emplace(v, VarBoundState::VAR_BOUND_AT_LOWER);
            } else {
                assert(false);
            }
            break;
        case QS_COL_BSTAT_UPPER:
            assert(simplex.hasUBound(v));
            newBounds.emplace(v, VarBoundState::VAR_BOUND_AT_UPPER);
            break;
        default:
            assert(false);
    }
}

LinResult ExternalQsoptex::qsxStatusToSolStatus() const {
    switch (d_status) {
        case QS_LP_OPTIMAL:
        case QS_LP_FEASIBLE:
        case QS_LP_DELTA_OPTIMAL:
        case QS_LP_DELTA_FEASIBLE:
        case QS_LP_UNBOUNDED:
            return isStrictVarZero() ? LinResult::LinInfeasible : LinResult::LinFeasible;
        case QS_LP_INFEASIBLE:
            return LinResult::LinInfeasible;
        case QS_LP_ITER_LIMIT:
        case QS_LP_TIME_LIMIT:
            return LinResult::LinExhausted;
        default:
            return LinResult::LinError;
    }
}

ExternalQsoptexEpsilon::ExternalQsoptexEpsilon(Simplex const & simplex, SMTConfig const & config)
    : ExternalQsoptex(simplex, config) {
    if (col_to_var.empty() || row_to_var.empty()) return;

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
        assert(tableau.isNonBasic(v));
        rat_class const & lb = get_bound<isEpsilon, true>(simplex, v, config.delta(), lb_temp);
        rat_class const & ub = get_bound<isEpsilon, false>(simplex, v, config.delta(), ub_temp);
        ensure_perturbed_bounds_do_not_overlap(*this, v, lb, ub, lb_temp, ub_temp);
        assert(lb <= ub);

        var_to_col.at(v.x) = numCols();
        rat_QSnew_col(d_qsx, rat_zeroLpNum, lb.data(), ub.data(),
#ifndef NDEBUG
                      (std::string("V") + std::to_string(v.x)).c_str()
#else
                      nullptr
#endif
        );
#ifndef NDEBUG
        rat_class temp;
        rat_QSget_bound(d_qsx, numCols() - 1, 'L', &*temp);
        assert(temp == lb);
        rat_QSget_bound(d_qsx, numCols() - 1, 'U', &*temp);
        assert(temp == ub);
#endif
    }

    // Construct the rows of the LP by parsing the polynomial constraints together
    // with the row bounds on the basic variables
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

        rat_class const & lb = get_bound<isEpsilon, true>(simplex, v, config.delta(), lb_temp);
        rat_class const & ub = get_bound<isEpsilon, false>(simplex, v, config.delta(), ub_temp);
        ensure_perturbed_bounds_do_not_overlap(*this, v, lb, ub, lb_temp, ub_temp);
        assert(lb <= ub);

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
void ExternalQsoptexEpsilon::solveImpl(unsigned int & precision, rat_class & delta) {
    [[maybe_unused]] int const res =
        QSdelta_solver_rat(d_qsx, *delta, reinterpret_cast<rat_t *>(d_x.data()), reinterpret_cast<rat_t *>(d_y.data()),
                           &d_basis, &precision, PRIMAL_SIMPLEX, &d_status, nullptr, nullptr);
    assert(res == 0);
}
Solution ExternalQsoptexEpsilon::extractSolution() {
    Solution sol(qsxStatusToSolStatus(), numCols(), numRows());
    // If we do not have a valid solution, just return the error back to the soplex.
    if (sol.linResult != LinResult::LinFeasible && sol.linResult != LinResult::LinInfeasible) return sol;
    assert(d_status == QS_LP_OPTIMAL || d_status == QS_LP_DELTA_OPTIMAL || d_status == QS_LP_FEASIBLE ||
           d_status == QS_LP_DELTA_FEASIBLE || d_status == QS_LP_UNBOUNDED || d_status == QS_LP_INFEASIBLE);

    for (int colIdx = 0; colIdx < numCols(); colIdx++) {
        extractVarValue<VariableType::COL>(colIdx, sol);
    }
    for (int rowIdx = 0; rowIdx < numRows(); rowIdx++) {
        extractVarValue<VariableType::ROW>(rowIdx, sol);
    }

    return sol;
}
ExternalQsoptexStrict::ExternalQsoptexStrict(Simplex const & simplex, SMTConfig const & config)
    : ExternalQsoptex(simplex, config) {
    if (col_to_var.empty() || row_to_var.empty()) return;

    Tableau const & tableau = simplex.getTableau();
    int const strict_var_col = static_cast<int>(col_to_var.size());

    // Construct the columns of the LP by assigning upper/lower bounds to each
    // variable
    std::uint32_t num_additional_strict_rows = 0;
    std::uint32_t nnz_strict_rows = 0;
    for (LVRef const v : col_to_var) {
        if (hasStrictLb(v)) {
            num_additional_strict_rows++;
            nnz_strict_rows += 2;
        }
        if (hasStrictUb(v)) {
            num_additional_strict_rows++;
            nnz_strict_rows += 2;
        }
    }
    for (LVRef const v : row_to_var) {
        if (hasStrictLb(v) || hasStrictUb(v)) {
            num_additional_strict_rows++;
            nnz_strict_rows += 2 + tableau.getRowPoly(v).size();
        }
    }

    std::vector<LVRef> rowToArithVarSplit;
    std::vector<int> numNonZeroPerRow;
    std::vector<int> beginRowIdx;
    std::vector<int> colIdxs;
    std::vector<rat_class> values;

    rowToArithVarSplit.reserve(row_to_var.size() + num_additional_strict_rows);
    d_rhs.reserve(row_to_var.size() + num_additional_strict_rows);
    d_range.reserve(row_to_var.size() + num_additional_strict_rows);
    numNonZeroPerRow.reserve(row_to_var.size() + num_additional_strict_rows);
    beginRowIdx.reserve(row_to_var.size() + num_additional_strict_rows);
    colIdxs.reserve(nnz + nnz_strict_rows);
    values.reserve(nnz + nnz_strict_rows);
    d_sense.reserve(row_to_var.size() + num_additional_strict_rows);

    rat_class lb_temp, ub_temp;
    for (LVRef const v : col_to_var) {
        assert(tableau.isNonBasic(v));
        rat_class const & lb = get_bound<isEpsilon, true>(simplex, v, config.delta(), lb_temp);
        rat_class const & ub = get_bound<isEpsilon, false>(simplex, v, config.delta(), ub_temp);

        bool const isLbStrict = hasStrictLb(v);
        bool const isUbStrict = hasStrictUb(v);

        if (isLbStrict) {
            assert(lb.isFinite());
            numNonZeroPerRow.emplace_back(2);
            beginRowIdx.emplace_back(colIdxs.size());
            colIdxs.emplace_back(numCols()); // Col of the variable we are about to add
            colIdxs.emplace_back(strict_var_col);
            values.emplace_back(1);
            values.emplace_back(-1);
            d_sense.emplace_back('G');
            d_rhs.emplace_back(lb);
            d_range.emplace_back(0);
            rowToArithVarSplit.emplace_back(v);
        }

        if (isUbStrict) {
            assert(ub.isFinite());
            numNonZeroPerRow.emplace_back(2);
            beginRowIdx.emplace_back(colIdxs.size());
            colIdxs.emplace_back(numCols()); // Col of the variable we are about to add
            colIdxs.emplace_back(strict_var_col);
            values.emplace_back(1);
            values.emplace_back(1);
            d_sense.emplace_back('L');
            d_rhs.emplace_back(ub);
            d_range.emplace_back(0);
            rowToArithVarSplit.emplace_back(v);
        }

        assert(lb <= ub);
        var_to_col.at(v.x) = numCols();
        rat_QSnew_col(d_qsx, rat_zeroLpNum, isLbStrict ? rat_NINFTY : lb.data(), isUbStrict ? rat_INFTY : ub.data(),
#ifndef NDEBUG
                      (std::string("V") + std::to_string(v.x)).c_str()
#else
                      nullptr
#endif
        );

#ifndef NDEBUG
        rat_class temp;
        rat_QSget_bound(d_qsx, numCols() - 1, 'L', &*temp);
        assert(temp == lb || hasStrictLb(v));
        rat_QSget_bound(d_qsx, numCols() - 1, 'U', &*temp);
        assert(temp == ub || hasStrictUb(v));
#endif
    }

    // At the end, add the strict variable column,
    // which is used to encode the strict bounds
    rat_QSnew_col(d_qsx, minus_one.data(), rat_zeroLpNum, rat_oneLpNum, nullptr);

    // Construct the rows of the LP by parsing the polynomial constraints together
    // with the row bounds on the auxiliary variables
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

        rat_class const & lb = get_bound<isEpsilon, true>(simplex, v, config.delta(), lb_temp);
        rat_class const & ub = get_bound<isEpsilon, false>(simplex, v, config.delta(), ub_temp);

        if (hasStrictBound(v)) {
            bool const isLbStrict = hasStrictLb(v);

            if (simplex.hasLBound(v)) {
                if (isLbStrict) {
                    numNonZeroPerRow.back()++;
                    colIdxs.emplace_back(strict_var_col);
                    values.emplace_back(-1);
                }
                assert(lb.isFinite());
                d_rhs.emplace_back(lb);
                d_range.emplace_back(0);
                d_sense.emplace_back('G');
                rowToArithVarSplit.emplace_back(v);
            }

            // Case II: we are dealing with a row with at least one bound. If both
            // bounds are set, we add the row twice, once for the lower bound and once
            // for the upper bound.
            if (hasBothBounds(v)) {
                // Insert all cols (except strict col, if present)
                beginRowIdx.emplace_back(colIdxs.size());
                colIdxs.insert(colIdxs.end(), colIdxs.end() - numNonZeroPerRow.back(),
                               isLbStrict ? colIdxs.end() - 1 : colIdxs.end());
                // Insert all values (except strict col, if present)
                values.insert(values.end(), values.end() - numNonZeroPerRow.back(),
                              isLbStrict ? values.end() - 1 : values.end());
                numNonZeroPerRow.emplace_back(isLbStrict ? numNonZeroPerRow.back() - 1 : numNonZeroPerRow.back());
            }

            if (simplex.hasUBound(v)) {
                if (hasStrictUb(v)) {
                    numNonZeroPerRow.back()++;
                    colIdxs.emplace_back(strict_var_col);
                    values.emplace_back(1);
                }
                assert(ub.isFinite());
                d_rhs.emplace_back(ub);
                d_range.emplace_back(0);
                d_sense.emplace_back('L');
                rowToArithVarSplit.emplace_back(v);
            }
        } else {
            assert(!hasStrictBound(v));
            rowToArithVarSplit.emplace_back(v);
            if (hasBothBounds(v)) // Bounded
            {
                assert(lb.isFinite() && ub.isFinite());
                d_rhs.emplace_back(lb);
                d_range.emplace_back(ub - lb);
                d_sense.emplace_back('R');
            } else if (simplex.hasLBound(v)) // Only lower bound
            {
                assert(lb.isFinite());
                d_rhs.emplace_back(lb);
                d_range.emplace_back(0);
                d_sense.emplace_back('G');
            } else if (simplex.hasUBound(v)) // Only upper bound
            {
                assert(ub.isFinite());
                d_rhs.emplace_back(ub);
                d_range.emplace_back(0);
                d_sense.emplace_back('L');
            } else // Unbounded
            {
                d_rhs.emplace_back(ninfinity);
                d_range.emplace_back(0);
                d_sense.emplace_back('G');
            }
        }
    }

    static_assert(sizeof(rat_class) == sizeof(rat_t), "mpq_class layout assumption broken");
    rat_QSadd_ranged_rows(d_qsx, static_cast<int>(numNonZeroPerRow.size()), numNonZeroPerRow.data(), beginRowIdx.data(),
                          colIdxs.data(), reinterpret_cast<rat_t const *>(values.data()),
                          reinterpret_cast<rat_t const *>(d_rhs.data()), d_sense.data(),
                          reinterpret_cast<rat_t const *>(d_range.data()), nullptr);

    row_to_var = std::move(rowToArithVarSplit);
    d_y.resize(row_to_var.size());
    assert(static_cast<int>(d_y.size()) == numRows());
}

void ExternalQsoptexStrict::solveImpl(unsigned int & precision, rat_class & delta) {
    rat_class obj_lb, obj_ub;
    [[maybe_unused]] int const res = QSdelta_full_solver_rat(
        d_qsx, *delta, reinterpret_cast<rat_t *>(d_x.data()), reinterpret_cast<rat_t *>(d_y.data()), *obj_lb, *obj_ub,
        &d_basis, PRIMAL_SIMPLEX, &d_status, &precision, nullptr, nullptr);
    assert(res == 0);
}
Solution ExternalQsoptexStrict::extractSolution() {
    Solution sol(qsxStatusToSolStatus(), numCols(), numRows());
    // If we do not have a valid solution, just return the error back to the soplex.
    if (sol.linResult != LinResult::LinFeasible && sol.linResult != LinResult::LinInfeasible) return sol;
    assert(d_status == QS_LP_OPTIMAL || d_status == QS_LP_DELTA_OPTIMAL || d_status == QS_LP_FEASIBLE ||
           d_status == QS_LP_DELTA_FEASIBLE || d_status == QS_LP_UNBOUNDED || d_status == QS_LP_INFEASIBLE);

    bool isStrictBasic = d_basis.cstat[numCols() - 1] == QS_COL_BSTAT_BASIC;

#ifndef NDEBUG
    std::cout << "Strict var value: " << d_x[numCols() - 1] << std::endl;
#endif

    // Set the t-var to 0 to simplify import. It is not needed if the problem is infeasible for some other reason
    if ((d_status == QS_LP_OPTIMAL || d_status == QS_LP_FEASIBLE || d_status == QS_LP_DELTA_OPTIMAL) && isStrictBasic &&
        d_x[numCols() - 1].isNonZero()) {
        int res = rat_QSchange_bound(d_qsx, numCols() - 1, 'U', rat_zeroLpNum);
        assert(res == 0);
        rat_class delta{config.delta() >= 0 ? config.delta() : 0.0};
        rat_class obj_lb, obj_ub;
        res = QSdelta_full_solver_rat(d_qsx, *delta, reinterpret_cast<rat_t *>(d_x.data()),
                                      reinterpret_cast<rat_t *>(d_y.data()), *obj_lb, *obj_ub, &d_basis, PRIMAL_SIMPLEX,
                                      &d_status, nullptr, nullptr, nullptr);
        assert(res == 0);
        assert(d_status == QS_LP_OPTIMAL || d_status == QS_LP_DELTA_OPTIMAL || d_status == QS_LP_FEASIBLE ||
               d_status == QS_LP_DELTA_FEASIBLE);
        isStrictBasic = d_basis.cstat[numCols() - 1] == QS_COL_BSTAT_BASIC;
    }

    // By this point, either the strict variable is basic and at 0,
    // or it is non-basic and we don't care about its value
    assert(!isStrictBasic || d_x[numCols() - 1].isZero());
    assert(d_status == QS_LP_OPTIMAL || d_status == QS_LP_DELTA_OPTIMAL || d_status == QS_LP_FEASIBLE ||
           d_status == QS_LP_DELTA_FEASIBLE || d_status == QS_LP_UNBOUNDED || d_status == QS_LP_INFEASIBLE);

    for (int rowIdx = 0; rowIdx < numRows(); rowIdx++) {
        extractVarValue<VariableType::ROW>(rowIdx, sol);
    }

    // This way, we skip the strict-var, since it is not included in col_to_var
    for (int colIdx = 0; colIdx < static_cast<int>(col_to_var.size()); colIdx++) {
        // We already encountered the nonbasic version of this variable, we can skip it
        if (sol.newNonBasis.contains(LVRef{col_to_var[colIdx]})) continue;
        extractVarValue<VariableType::COL>(colIdx, sol);
    }

    // Since we have to split range rows, we have to remove those that we
    // erroneously marked as basic and then discovered to be nonbasic
    // for (LVRef const v : sol.newNonBasis) {
    //     sol.newBasis.erase(v);
    // }
    // If the strict variable is basic, we need to add some other non-basic
    // variable to the basis to maintain the same number of basic variables
    return sol;
}

std::unique_ptr<ExternalLpSolver> make_external_qsoptex(Simplex const & simplex, SMTConfig const & config) {
    if (config.epsilon_mode()) return std::make_unique<ExternalQsoptexEpsilon>(simplex, config);
    return std::make_unique<ExternalQsoptexStrict>(simplex, config);
}

} // namespace opensmt::externallp

#else

namespace opensmt::externallp {
std::unique_ptr<ExternalLpSolver> make_external_qsoptex(Simplex const &, SMTConfig const &) {
    throw std::runtime_error("opensmt was not build with support for the qsoptex lp solver");
}
} // namespace opensmt::externallp

#endif
