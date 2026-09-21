#ifndef DM_SHIFT_H
#define DM_SHIFT_H

#include <complex>

namespace module_dm
{
/**
 * @brief map a real/complex type to the opposite one
 * ShiftRealComplex<double>::type = std::complex<double>
 * ShiftRealComplex<std::complex<double>>::type = double
 */
template<typename T> struct ShiftRealComplex
{
    using type = void;
};

template<>
struct ShiftRealComplex<double>
{
    using type = std::complex<double>;
};

template<>
struct ShiftRealComplex<std::complex<double>>
{
    using type = double;
};

} // namespace module_dm

#endif
