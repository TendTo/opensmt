#include "ExternalSoplex.h"

#ifdef OPENSSMT_ENABLE_SOPLEX

namespace opensmt::externallp {

using SoPlex = soplex::SoPlex;

namespace {
    class Rational : public soplex::Rational {
    public:
        using soplex::Rational::Rational;
        using soplex::Rational::operator=;
        using soplex::Rational::compare;

        static Rational const minus_one;
        static Rational const infinity;
        static Rational const ninfinity;

        explicit Rational(FastRational const & val) : Rational() { *this = val; }
        explicit Rational(std::int32_t const val) : Rational() { *this = val; }

        Rational & operator=(FastRational const & val) {
            if (auto const pair = val.tryGetNumDen(); pair.has_value()) {
                backend().set(pair.value().first, pair.value().second);
            } else {
                backend() = val.getMpq();
            }
            return *this;
        }
        Rational & operator=(Delta const & val) { return *this = val.R(); }

        std::strong_ordering operator<=>(Rational const & o) const {
            int const res = compare(o);
            if (res == 0) return std::strong_ordering::equal;
            return res < 0 ? std::strong_ordering::less : std::strong_ordering::greater;
        }

        operator Delta() const { return Delta{static_cast<FastRational>(*this)}; }
        operator FastRational() const {
            return backend().is_ratgmp() ? FastRational(reinterpret_cast<mpq_class &>(*backend().get_gmp()))
                                         : FastRational{backend().num(), backend().den()};
        }

        Rational & add(std::int32_t const num, std::int32_t const den) {
            *this += soplex::Rational(num, den);
            return *this;
        }
        Rational & sub(std::int32_t const num, std::int32_t const den) { return add(-num, den); }
    };

    Rational const Rational::infinity{std::numeric_limits<double>::infinity()};
    Rational const Rational::ninfinity{-std::numeric_limits<double>::infinity()};
    Rational const Rational::minus_one{-1};

    std::strong_ordering operator<=>(soplex::Rational const & r, FastRational const & o) {
        int res;
        if (auto const pair = o.tryGetNumDen(); pair.has_value()) {
            res = r.backend().compare(pair.value().first, pair.value().second);
        } else {
            res = r.backend().compare(o.getMpq());
        }
        if (res == 0) return std::strong_ordering::equal;
        return res < 0 ? std::strong_ordering::less : std::strong_ordering::greater;
    }

    FastRational to_fast_rational(soplex::Rational const & r) {
        return r.backend().is_rat32() ? FastRational(r.backend().num(), r.backend().den())
                                      : FastRational(reinterpret_cast<mpq_class &>(*r.backend().get_gmp()));
    }

    template<bool epsilon_mode, bool is_lower_bound>
    Rational const & get_bound(Simplex const & simplex, LVRef const v, double const delta, Rational & bound) {
        if (bool const is_finite = is_lower_bound ? simplex.hasLBound(v) : simplex.hasUBound(v); !is_finite)
            return is_lower_bound ? Rational::ninfinity : Rational::infinity;

        assert(is_lower_bound ? simplex.hasLBound(v) : simplex.hasUBound(v));
        Delta const & delta_bound = is_lower_bound ? simplex.Lb(v) : simplex.Ub(v);
        bool const is_strict = delta < 0 && delta_bound.hasDelta();

        bound = delta_bound.R();
        if (epsilon_mode && is_strict) {
            // Adjust the bounds with epsilon
            if (is_lower_bound) {
                bound.add(ExternalLpSolver::epsilon_num, ExternalLpSolver::epsilon_den);
            } else {
                bound.sub(ExternalLpSolver::epsilon_num, ExternalLpSolver::epsilon_den);
            }
        }
        return bound;
    }

    void ensure_perturbed_bounds_do_not_overlap(ExternalLpSolver const & solver, LVRef v, Rational const & lb,
                                                Rational const & ub, Rational & lb_temp, Rational & ub_temp) {
        if (solver.hasStrictBound(v) && lb > ub) {
            assert(lb.backend().is_finite() && ub.backend().is_finite());
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

ExternalSoplex::ExternalSoplex(Simplex const & simplex, SMTConfig const & config)
    : ExternalLpSolver(simplex, config),
      d_spx() {
    d_spx.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
    d_spx.setIntParam(SoPlex::SIMPLIFIER, SoPlex::SIMPLIFIER_OFF);
    d_spx.setIntParam(SoPlex::ALGORITHM, SoPlex::ALGORITHM_PRIMAL);
    d_spx.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
    d_spx.setIntParam(SoPlex::READMODE, SoPlex::READMODE_RATIONAL);
    d_spx.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
    d_spx.setIntParam(SoPlex::CHECKMODE, SoPlex::CHECKMODE_RATIONAL);
    d_spx.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
    d_spx.setIntParam(SoPlex::PRICER, SoPlex::PRICER_AUTO);
    d_spx.setIntParam(SoPlex::ITERLIMIT, lp_pivot_threshold);
    // If delta, allow for feasibility tolerance
    d_spx.setRealParam(SoPlex::FEASTOL, std::max(config.delta(), 0.0));
    d_spx.setRealParam(SoPlex::OPTTOL, 0.0);

    if (config.verbosity() >= 2) d_spx.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_DEBUG);
}

Solution ExternalSoplex::solve() {
    if (col_to_var.empty() || row_to_var.empty()) return {LinResult::LinEmtpy};

    SpxStatus res = SpxStatus::UNKNOWN;
    try {
        res = d_spx.optimize();
    } catch (soplex::SPxException const &) { return {LinResult::LinError}; }

    stats.iterative_refinements[d_spx.numRefinements()]++;
    std::size_t precision = d_spx.numPrecisionBoosts() == 0 ? sizeof(double) * 8 : 167;
    for (int i = 1; i < d_spx.numPrecisionBoosts(); i++) {
        precision = static_cast<std::size_t>(precision * d_spx.realParam(SoPlex::PRECISION_BOOSTING_FACTOR));
    }
    stats.precisions[precision]++;
    if (config.delta() >= 0) {
        soplex::Rational maxBoundViolation, sumBoundViolation;
        soplex::Rational maxRowViolation, sumRowViolation;
        d_spx.getBoundViolationRational(maxBoundViolation, sumBoundViolation);
        d_spx.getRowViolationRational(maxRowViolation, sumRowViolation);
        double const delta = std::max(maxBoundViolation.convert_to<double>(), maxRowViolation.convert_to<double>());
        stats.max_delta = std::max(stats.max_delta, delta);
        if (delta > 0) stats.delta_calls++;
    }

    [[maybe_unused]] bool status;
    switch (res) {
        case SpxStatus::OPTIMAL:
        case SpxStatus::UNBOUNDED:
            assert(d_spx.hasSol());
            d_x.reDim(numCols());
            status = d_spx.getPrimalRational(d_x);
            assert(status);
            // Check the value of the last column (strict variable)
            d_sol_status = isStrictVarZero() ? LinResult::LinInfeasible : LinResult::LinFeasible;
            return extractSolution();
        case SpxStatus::INFEASIBLE:
            assert(d_spx.hasDualFarkas());
            d_y.reDim(numRows());
            status = d_spx.getDualRational(d_y);
            assert(status);
            d_sol_status = LinResult::LinInfeasible;
            return extractSolution();
        case SpxStatus::ABORT_ITER:
        case SpxStatus::ABORT_TIME:
        case SpxStatus::ABORT_CYCLING:
            return {LinResult::LinExhausted};
        default:
            return {LinResult::LinError};
    }
}

template<ExternalSoplex::VariableType VarType>
void ExternalSoplex::extractVarValue(int const idx, Solution & sol, soplex::VectorRational const * const values) {
    auto & newNonBasis = sol.newNonBasis;
    auto & newValues = sol.newValues;
    auto & newBounds = sol.newBounds;

    LVRef v = LVRef::Undef;
    VarStatus varStatus = VarStatus::UNDEFINED;

    if constexpr (VarType == VariableType::COL) {
        v = col_to_var.at(idx);
        varStatus = d_spx.basisColStatus(idx);
    }
    if constexpr (VarType == VariableType::ROW) {
        v = row_to_var.at(idx);
        varStatus = d_spx.basisRowStatus(idx);
    }
    assert(v != LVRef::Undef);

    if (varStatus != VarStatus::BASIC) newNonBasis.emplace(v);

    soplex::Rational value_tmp;
    soplex::Rational const * value = nullptr;
    switch (varStatus) {
        // If we are dealing with a basic variable, we necessarily need to get its
        // value from the solved problem.
        case VarStatus::BASIC:
        case VarStatus::UNDEFINED:
            if (values != nullptr) {
                value = &(*values)[idx];
            } else if (VarType == VariableType::ROW) {
                d_spx.getRowActivityRational(idx, value_tmp);
                value = &value_tmp;
            } else if (VarType == VariableType::COL) {
                d_spx.getColActivityRational(idx, value_tmp);
                value = &value_tmp;
            }
            assert(value != nullptr);
            if (simplex.hasLBound(v) && *value <= simplex.Lb(v).R()) {
                newBounds.emplace(v, VarBoundState::VAR_BOUND_AT_LOWER);
            } else if (simplex.hasUBound(v) && *value >= simplex.Ub(v).R()) {
                newBounds.emplace(v, VarBoundState::VAR_BOUND_AT_UPPER);
            } else {
                newValues.emplace(v, Delta(to_fast_rational(*value)));
            }
            break;
        // For non-basic variables we know the value is at a bound
        // and we can use the d_vars/proper LP bounds directly
        case VarStatus::FIXED:
            assert(simplex.hasLBound(v));
            newBounds.insert_or_assign(v, VarBoundState::VAR_BOUND_AT_LOWER);
            break;
        case VarStatus::ON_LOWER:
            assert(simplex.hasLBound(v));
            newBounds.insert_or_assign(v, VarBoundState::VAR_BOUND_AT_LOWER);
            break;
        case VarStatus::ON_UPPER:
            assert(simplex.hasUBound(v));
            newBounds.insert_or_assign(v, VarBoundState::VAR_BOUND_AT_UPPER);
            break;
        case VarStatus::ZERO:
            newValues.emplace(v, Delta());
            break;
        default:
            assert(false);
    }
}

ExternalSoplexEpsilon::ExternalSoplexEpsilon(Simplex const & simplex, SMTConfig const & config)
    : ExternalSoplex(simplex, config) {
    if (col_to_var.empty() || row_to_var.empty()) return;

    Tableau const & tableau = simplex.getTableau();

    // The number of cols must accommodate for the non-aux variables as well as
    // the additional strict variable t
    soplex::LPRowSetRational rows(static_cast<int>(row_to_var.size()));
    soplex::LPColSetRational cols(static_cast<int>(col_to_var.size()));

    Rational lb_temp, ub_temp;
    for (LVRef const v : col_to_var) {
        assert(tableau.isNonBasic(v));
        Rational const & lb = get_bound<isEpsilon, true>(simplex, v, config.delta(), lb_temp);
        Rational const & ub = get_bound<isEpsilon, false>(simplex, v, config.delta(), ub_temp);
        ensure_perturbed_bounds_do_not_overlap(*this, v, lb, ub, lb_temp, ub_temp);
        assert(lb <= ub);

        var_to_col.at(v.x) = cols.num();
        cols.add({0.0, soplex::DSVectorRational(), ub, lb});
    }

    // Construct the rows of the LP by parsing the polynomial constraints together
    // with the row bounds on the basic variables
    for (LVRef const v : row_to_var) {
        assert(tableau.isBasic(v));

        PolynomialT<LVRef> const & p = tableau.getRowPoly(v);
        soplex::DSVectorRational vec(static_cast<int>(p.size()));
        for (auto const & term : p) {
            int const colIdx = var_to_col[term.var.x];
            assert(colIdx >= 0);
            vec.add(colIdx, Rational{term.coeff});
        }

        Rational const & lb = get_bound<isEpsilon, true>(simplex, v, config.delta(), lb_temp);
        Rational const & ub = get_bound<isEpsilon, false>(simplex, v, config.delta(), ub_temp);
        ensure_perturbed_bounds_do_not_overlap(*this, v, lb, ub, lb_temp, ub_temp);
        assert(lb <= ub);

        rows.add({lb, vec, ub});
    }

    // Add both columns and rows to the LP
    d_spx.addColsRational(cols);
    d_spx.addRowsRational(rows);
}

Solution ExternalSoplexEpsilon::extractSolution() {
    Solution sol(d_sol_status, numCols(), numRows());
    // If we do not have a valid solution, just return the error back to the soplex.
    if (sol.linResult != LinResult::LinFeasible && sol.linResult != LinResult::LinInfeasible) return sol;
    assert(d_spx.status() == SpxStatus::INFEASIBLE || d_spx.status() == SpxStatus::OPTIMAL ||
           d_spx.status() == SpxStatus::UNBOUNDED);

    for (int colIdx = 0; colIdx < numCols(); colIdx++) {
        extractVarValue<VariableType::COL>(colIdx, sol, d_spx.status() != SpxStatus::INFEASIBLE ? &d_x : nullptr);
    }
    for (int rowIdx = 0; rowIdx < numRows(); rowIdx++) {
        extractVarValue<VariableType::ROW>(rowIdx, sol, d_spx.status() == SpxStatus::INFEASIBLE ? &d_y : nullptr);
    }

    return sol;
}
ExternalSoplexStrict::ExternalSoplexStrict(Simplex const & simplex, SMTConfig const & config)
    : ExternalSoplex(simplex, config) {
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

    soplex::LPRowSetRational rows(static_cast<int>(row_to_var.size() + num_additional_strict_rows));
    soplex::LPColSetRational cols(static_cast<int>(col_to_var.size()));
    std::vector<LVRef> rowToArithVarSplit;
    rowToArithVarSplit.reserve(row_to_var.size() + num_additional_strict_rows);

    Rational lb_temp, ub_temp;
    for (LVRef const v : col_to_var) {
        assert(tableau.isNonBasic(v));
        Rational const & lb = get_bound<isEpsilon, true>(simplex, v, config.delta(), lb_temp);
        Rational const & ub = get_bound<isEpsilon, false>(simplex, v, config.delta(), ub_temp);

        bool const isLbStrict = hasStrictLb(v);
        bool const isUbStrict = hasStrictUb(v);

        if (isLbStrict) {
            assert(lb.backend().is_finite());
            soplex::DSVectorRational vec(2);
            vec.add(cols.num(), 1);
            vec.add(strict_var_col, -1);
            rows.add({lb, vec, soplex::infinity});
            rowToArithVarSplit.emplace_back(v);
        }

        if (isUbStrict) {
            assert(ub.backend().is_finite());
            soplex::DSVectorRational vec(2);
            vec.add(cols.num(), 1);
            vec.add(strict_var_col, 1);
            rows.add({-soplex::infinity, vec, ub});
            rowToArithVarSplit.emplace_back(v);
        }

        assert(lb <= ub);
        var_to_col.at(v.x) = cols.num();
        cols.add({0.0, soplex::DSVectorRational(), isUbStrict ? Rational::infinity : ub,
                  isLbStrict ? Rational::ninfinity : lb});
    }

    // At the end, add the strict variable column,
    // which is used to encode the strict bounds
    cols.add({-1, soplex::DSVectorRational(), 1, 0});

    // Construct the rows of the LP by parsing the polynomial constraints together
    // with the row bounds on the auxiliary variables
    for (LVRef const v : row_to_var) {
        assert(tableau.isBasic(v));

        PolynomialT<LVRef> const & p = tableau.getRowPoly(v);
        soplex::DSVectorRational vec(static_cast<int>(p.size()) + 1);
        for (auto const & term : p) {
            int const colIdx = var_to_col[term.var.x];
            assert(colIdx >= 0);
            vec.add(colIdx, Rational{term.coeff});
        }

        Rational const & lb = get_bound<isEpsilon, true>(simplex, v, config.delta(), lb_temp);
        Rational const & ub = get_bound<isEpsilon, false>(simplex, v, config.delta(), ub_temp);

        if (hasStrictBound(v)) {
            if (simplex.hasLBound(v)) {
                if (hasStrictLb(v)) vec.add(strict_var_col, -1);
                assert(lb.backend().is_finite());
                rows.add({lb, vec, Rational::infinity});
                rowToArithVarSplit.emplace_back(v);
            }

            // Case II: we are dealing with a row with at least one bound. If both
            // bounds are set, we add the row twice, once for the lower bound and once
            // for the upper bound.
            if (simplex.hasUBound(v)) {
                // Ensure that the strict variable is present only once and with the
                // correct coefficient in the row vector
                if (int const idx = vec.pos(strict_var_col); idx > -1) vec.remove(idx);
                if (hasStrictUb(v)) vec.add(strict_var_col, 1);
                assert(ub.backend().is_finite());
                rows.add({Rational::ninfinity, vec, ub});
                rowToArithVarSplit.emplace_back(v);
            }
        } else {
            assert(!hasStrictBound(v));
            rowToArithVarSplit.emplace_back(v);
            rows.add({lb, vec, ub});
        }
    }

    // Add both columns and rows to the LP
    d_spx.addColsRational(cols);
    d_spx.addRowsRational(rows);
    row_to_var = std::move(rowToArithVarSplit);
}

Solution ExternalSoplexStrict::extractSolution() {
    Solution sol(d_sol_status, numCols(), numRows());
    // If we do not have a valid solution, just return the error back to the soplex.
    if (sol.linResult != LinResult::LinFeasible && sol.linResult != LinResult::LinInfeasible) return sol;
    assert(d_spx.status() == SpxStatus::INFEASIBLE || d_spx.status() == SpxStatus::OPTIMAL);

    bool const isStrictNonBasic = d_spx.basisColStatus(d_spx.numColsRational() - 1) != VarStatus::BASIC;

    // Set the t-var to 0 to simplify import when the t-var is non-basic.
    // It is not needed if the problem is infeasible for some other reason
    if (d_spx.status() == SpxStatus::OPTIMAL && isStrictNonBasic && !d_x[numCols() - 1].is_zero()) {
        assert(d_spx.hasSol());
        assert(d_x.dim() == d_spx.numColsRational());
        d_spx.changeBoundsRational(d_spx.numColsRational() - 1, soplex::Rational(), soplex::Rational());
        d_spx.optimize();
        assert(d_spx.status() == SpxStatus::OPTIMAL);

        assert(d_x.dim() == d_spx.numColsRational());
        [[maybe_unused]] bool const getPrimalSuccess = d_spx.getPrimalRational(d_x);
        assert(getPrimalSuccess);
    }

#ifndef NDEBUG
    if (d_spx.status() == SpxStatus::OPTIMAL) std::cout << "Strict var value: " << d_x[numCols() - 1] << std::endl;
#endif

    // By this point, either the strict variable is non-basic and at 0,
    // or it is basic and we don't care about its value
    assert(d_spx.status() == SpxStatus::INFEASIBLE || d_spx.status() == SpxStatus::OPTIMAL);
    assert(d_spx.status() == SpxStatus::INFEASIBLE || !isStrictNonBasic || d_x[numCols() - 1].is_zero());

    for (int rowIdx = 0; rowIdx < numRows(); rowIdx++) {
        extractVarValue<VariableType::ROW>(rowIdx, sol, d_spx.status() == SpxStatus::INFEASIBLE ? &d_y : nullptr);
    }

    // This way, we skip the strict-var, since it is not included in col_to_var
    for (int colIdx = 0; colIdx < static_cast<int>(col_to_var.size()); colIdx++) {
        // We already encountered the nonbasic version of this variable, we can skip it
        if (sol.newNonBasis.contains(LVRef{col_to_var[colIdx]})) continue;
        extractVarValue<VariableType::COL>(colIdx, sol, d_spx.status() != SpxStatus::INFEASIBLE ? &d_x : nullptr);
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

std::unique_ptr<ExternalLpSolver> make_external_soplex(Simplex const & simplex, SMTConfig const & config) {
    if (config.epsilon_mode()) return std::make_unique<ExternalSoplexEpsilon>(simplex, config);
    return std::make_unique<ExternalSoplexStrict>(simplex, config);
}
} // namespace opensmt::externallp

#else

namespace opensmt::externallp {
std::unique_ptr<ExternalLpSolver> make_external_soplex(Simplex const &, SMTConfig const &) {
    throw std::runtime_error("opensmt was not build with support for the soplex lp solver");
}
} // namespace opensmt::externallp

#endif
