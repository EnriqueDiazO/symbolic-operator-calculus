"""Safe substitution of free SymPy variables for MVP integral expressions."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import sympy as sp


class SafeSubstitutionError(ValueError):
    """Raised when a safe free-variable substitution cannot be guaranteed."""


DEFAULT_DUMMY_NAMES = ("y", "u", "v", "w", "z")


def collect_bound_symbols(expression: sp.Expr) -> set[sp.Symbol]:
    """Collect bound integration variables from all nested SymPy Integrals."""

    if isinstance(expression, sp.Integral):
        bound = set(expression.variables)
        bound.update(collect_bound_symbols(expression.function))
        for limit in expression.limits:
            for part in limit[1:]:
                bound.update(collect_bound_symbols(part))
        return bound

    bound: set[sp.Symbol] = set()
    for arg in expression.args:
        bound.update(collect_bound_symbols(arg))
    return bound


def fresh_symbol(
    avoid: Iterable[sp.Symbol],
    *,
    preferred_names: Sequence[str] = DEFAULT_DUMMY_NAMES,
) -> sp.Symbol:
    """Return a deterministic fresh Symbol avoiding existing symbol names."""

    avoid_names = {symbol.name for symbol in avoid}
    for name in preferred_names:
        if name not in avoid_names:
            return sp.Symbol(name)

    index = 1
    while True:
        for name in preferred_names:
            candidate = f"{name}{index}"
            if candidate not in avoid_names:
                return sp.Symbol(candidate)
        index += 1


def substitute_free_variable(
    expression: sp.Expr,
    old_variable: sp.Symbol,
    new_expression: sp.Expr,
) -> sp.Expr:
    """Substitute only free occurrences of ``old_variable``.

    Bound variables of ``Integral`` are preserved. If the replacement contains
    a symbol that would be captured by an integral, that bound variable is
    alpha-renamed first using a deterministic fresh symbol.
    """

    if not isinstance(expression, sp.Expr):
        raise SafeSubstitutionError("expression must be a SymPy expression.")
    if not isinstance(old_variable, sp.Symbol):
        raise SafeSubstitutionError("old_variable must be a SymPy Symbol.")
    if not isinstance(new_expression, sp.Expr):
        raise SafeSubstitutionError("new_expression must be a SymPy expression.")
    return _substitute_free(expression, old_variable, new_expression)


def _substitute_free(
    expression: sp.Expr,
    old_variable: sp.Symbol,
    new_expression: sp.Expr,
) -> sp.Expr:
    if expression == old_variable:
        return new_expression
    if not expression.args:
        return expression

    if isinstance(expression, sp.Integral):
        return _substitute_free_in_integral(expression, old_variable, new_expression)

    new_args = tuple(
        _substitute_free(arg, old_variable, new_expression)
        for arg in expression.args
    )
    if new_args == expression.args:
        return expression
    return expression.func(*new_args)


def _substitute_free_in_integral(
    integral: sp.Integral,
    old_variable: sp.Symbol,
    new_expression: sp.Expr,
) -> sp.Integral:
    bound_variables = set(integral.variables)
    new_limits = _substitute_free_in_limits(
        integral.limits,
        old_variable,
        new_expression,
    )

    if old_variable in bound_variables:
        return sp.Integral(integral.function, *new_limits)

    prepared = integral
    conflicting = bound_variables.intersection(new_expression.free_symbols)
    if conflicting:
        avoid = (
            integral.free_symbols
            | collect_bound_symbols(integral)
            | new_expression.free_symbols
            | new_expression.atoms(sp.Symbol)
            | {old_variable}
        )
        rename_map: dict[sp.Symbol, sp.Symbol] = {}
        for symbol in integral.variables:
            if symbol in conflicting:
                fresh = fresh_symbol(avoid)
                rename_map[symbol] = fresh
                avoid.add(fresh)
        prepared = alpha_rename_integral(integral, rename_map)
        new_limits = _substitute_free_in_limits(
            prepared.limits,
            old_variable,
            new_expression,
        )

    new_function = _substitute_free(prepared.function, old_variable, new_expression)
    return sp.Integral(new_function, *new_limits)


def _substitute_free_in_limits(
    limits: tuple[tuple[sp.Expr, ...], ...],
    old_variable: sp.Symbol,
    new_expression: sp.Expr,
) -> tuple[tuple[sp.Expr, ...], ...]:
    new_limits: list[tuple[sp.Expr, ...]] = []
    for limit in limits:
        if len(limit) == 1:
            new_limits.append(limit)
            continue
        variable = limit[0]
        new_parts = tuple(
            _substitute_free(part, old_variable, new_expression)
            for part in limit[1:]
        )
        new_limits.append((variable, *new_parts))
    return tuple(new_limits)


def alpha_rename_integral(
    integral: sp.Integral,
    rename_map: dict[sp.Symbol, sp.Symbol],
) -> sp.Integral:
    """Alpha-rename variables bound by this Integral only."""

    if not rename_map:
        return integral

    renamed_function = integral.function
    for old, new in rename_map.items():
        renamed_function = _rename_bound_occurrences(renamed_function, old, new)

    renamed_limits: list[tuple[sp.Expr, ...]] = []
    for limit in integral.limits:
        variable = limit[0]
        new_variable = rename_map.get(variable, variable)
        renamed_parts: list[sp.Expr] = []
        for part in limit[1:]:
            renamed_part = part
            for old, new in rename_map.items():
                renamed_part = _rename_bound_occurrences(renamed_part, old, new)
            renamed_parts.append(renamed_part)
        renamed_limits.append((new_variable, *renamed_parts))
    return sp.Integral(renamed_function, *renamed_limits)


def _rename_bound_occurrences(
    expression: sp.Expr,
    old: sp.Symbol,
    new: sp.Symbol,
) -> sp.Expr:
    if expression == old:
        return new
    if not expression.args:
        return expression

    if isinstance(expression, sp.Integral):
        inner_bound = set(expression.variables)
        if old in inner_bound:
            renamed_limits = []
            for limit in expression.limits:
                variable = limit[0]
                renamed_parts = tuple(
                    _rename_bound_occurrences(part, old, new)
                    for part in limit[1:]
                )
                renamed_limits.append((variable, *renamed_parts))
            return sp.Integral(expression.function, *renamed_limits)
        renamed_function = _rename_bound_occurrences(expression.function, old, new)
        renamed_limits = []
        for limit in expression.limits:
            variable = limit[0]
            renamed_parts = tuple(
                _rename_bound_occurrences(part, old, new)
                for part in limit[1:]
            )
            renamed_limits.append((variable, *renamed_parts))
        return sp.Integral(renamed_function, *renamed_limits)

    new_args = tuple(
        _rename_bound_occurrences(arg, old, new)
        for arg in expression.args
    )
    if new_args == expression.args:
        return expression
    return expression.func(*new_args)
