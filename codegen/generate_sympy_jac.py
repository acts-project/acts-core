import sympy as sym
from sympy import MatrixSymbol

from sympy_common import name_expr, find_by_name, cxx_printer, my_expression_print


step_path_derivatives = (
    MatrixSymbol("step_path_derivatives", 8, 1).as_explicit().as_mutable()
)
step_path_derivatives[7, 0] = 0  # qop

surface_path_derivatives = (
    MatrixSymbol("surface_path_derivatives", 1, 8).as_explicit().as_mutable()
)
surface_path_derivatives[0, 3] = 0
surface_path_derivatives[0, 7] = 0

J_bf = MatrixSymbol("J_bf", 8, 6).as_explicit().as_mutable()
tmp = sym.zeros(8, 6)
tmp[0:3, 0:2] = J_bf[0:3, 0:2]
tmp[0:3, 2:4] = J_bf[0:3, 2:4]
tmp[4:7, 2:4] = J_bf[4:7, 2:4]  # line surface
tmp[3, 5] = 1
tmp[7, 4] = 1
J_bf = tmp

J_t = MatrixSymbol("J_t", 8, 8).as_explicit().as_mutable()
tmp = sym.eye(8)
tmp[0:3, 4:8] = J_t[0:3, 4:8]
tmp[3, 7] = J_t[3, 7]
tmp[4:7, 4:8] = J_t[4:7, 4:8]
J_t = tmp

J_fb = MatrixSymbol("J_fb", 6, 8).as_explicit().as_mutable()
tmp = sym.zeros(6, 8)
tmp[0:2, 0:3] = J_fb[0:2, 0:3]
tmp[2:4, 4:7] = J_fb[2:4, 4:7]
tmp[5, 3] = 1
tmp[4, 7] = 1
J_fb = tmp


def full_transport_jacobian_generic():
    J_full = name_expr(
        "J_full",
        J_fb
        * (sym.eye(8) + step_path_derivatives * surface_path_derivatives)
        * J_t
        * J_bf,
    )

    return [J_full]


def full_transport_jacobian_curvilinear(direction):
    surface_path_derivatives = (
        MatrixSymbol("surface_path_derivatives", 1, 8).as_explicit().as_mutable()
    )
    surface_path_derivatives[0, 0:3] = -direction.as_explicit().transpose()
    surface_path_derivatives[0, 3:8] = sym.zeros(1, 5)

    J_full = name_expr(
        "J_full",
        J_fb
        * (sym.eye(8) + step_path_derivatives * surface_path_derivatives)
        * J_t
        * J_bf,
    )

    return [J_full]


def my_full_transport_jacobian_generic_function_print(name_exprs, run_cse=True):
    printer = cxx_printer
    outputs = [find_by_name(name_exprs, name)[0] for name in ["J_full"]]

    lines = []

    head = "template <typename T> void boundToBoundTransportJacobianImpl(const T* J_fb, const T* J_t, const T* J_bf, const T* step_path_derivatives, const T* surface_path_derivatives, T* J_full) {"
    lines.append(head)

    code = my_expression_print(
        printer,
        name_exprs,
        outputs,
        run_cse=run_cse,
    )
    lines.extend([f"  {l}" for l in code.split("\n")])

    lines.append("}")

    return "\n".join(lines)


def my_full_transport_jacobian_curvilinear_function_print(name_exprs, run_cse=True):
    printer = cxx_printer
    outputs = [find_by_name(name_exprs, name)[0] for name in ["J_full"]]

    lines = []

    head = "template <typename T> void boundToCurvilinearTransportJacobianImpl(const T* J_fb, const T* J_t, const T* J_bf, const T* step_path_derivatives, const T* dir, T* J_full) {"
    lines.append(head)

    code = my_expression_print(
        printer,
        name_exprs,
        outputs,
        run_cse=run_cse,
    )
    lines.extend([f"  {l}" for l in code.split("\n")])

    lines.append("}")

    return "\n".join(lines)


print(
    """// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Note: This file is generated by generate_sympy_jac.py
//       Do not modify it manually.

#pragma once

#include <cmath>
"""
)

all_name_exprs = full_transport_jacobian_generic()
code = my_full_transport_jacobian_generic_function_print(
    all_name_exprs,
    run_cse=True,
)
print(code)
print()

all_name_exprs = full_transport_jacobian_curvilinear(MatrixSymbol("dir", 3, 1))
code = my_full_transport_jacobian_curvilinear_function_print(
    all_name_exprs,
    run_cse=True,
)
print(code)
