#include "ExternalQsoptex.h"

#include <utility>

#ifdef OPENSSMT_ENABLE_QSOPTEX

namespace opensmt::externallp {

namespace {

    constexpr int pivot_threshold = 1000000;

    rat_class minus_one;

    void QSXStart() {
        static bool started = false;
        if (started) return;
        started = true;
        QSexactStart();
        minus_one = -1;
    }

    void QSXFinish() {
        QSexactClear();
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
rat_class::operator FastRational() const {
    return rat_is_gmp(rat) ? FastRational(static_cast<mpq_class>(rat_get_mpq(rat)))
                           : FastRational{rat_get_num(rat), rat_get_den(rat)};
}

ExternalQsoptex::ExternalQsoptex(SMTConfig const & config) : ExternalLpSolver(config), d_status(), d_basis{0} {
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
