// SPDX-License-Identifier: MIT
// Copyright (c) 2019-2021 Thomas Vanderbruggen <th.vanderbruggen@gmail.com>

#ifndef SCICPP_CORE_UTILS
#define SCICPP_CORE_UTILS

#include "scicpp/core/meta.hpp"

#include <array>
#include <iterator>
#include <vector>

namespace scicpp::utils {

//---------------------------------------------------------------------------------
// set_array: Set an array of same size than a given array
//---------------------------------------------------------------------------------

template <typename OutputType, typename T, std::size_t N>
auto set_array(const std::array<T, N> & /* unused */) {
    return std::array<OutputType, N>{};
}

template <typename OutputType, typename T>
auto set_array(std::vector<T> v) {
    return std::vector<OutputType>(v.size());
}

template <class Array>
auto set_array(const Array &a) {
    return set_array<typename Array::value_type>(a);
}

//---------------------------------------------------------------------------------
// subvector
//---------------------------------------------------------------------------------

template <typename T,
          typename DiffTp = typename std::vector<T>::difference_type>
auto subvector(std::vector<T> &&v, std::size_t len, DiffTp offset = 0) {
    return std::vector<T>(
        std::make_move_iterator(v.begin() + offset),
        std::make_move_iterator(v.begin() + offset +
                                int(std::min(len, v.size()))));
}

template <typename Array, typename DiffTp = typename Array::difference_type>
auto subvector(const Array &v, std::size_t len, DiffTp offset = 0) {
    using T = typename Array::value_type;
    static_assert(meta::is_iterable_v<Array>);

    return std::vector<T>(
        std::make_move_iterator(v.begin() + offset),
        std::make_move_iterator(v.begin() + offset +
                                int(std::min(len, v.size()))));
}

//---------------------------------------------------------------------------------
// initialize zero
//---------------------------------------------------------------------------------

template <typename T>
constexpr auto set_zero() {
    if constexpr (meta::is_complex_v<T>) {
        using U = typename T::value_type;
        return T(U(0), U(0));
    } else {
        return T(0);
    }
}

} // namespace scicpp::utils

#endif // SCICPP_CORE_UTILS