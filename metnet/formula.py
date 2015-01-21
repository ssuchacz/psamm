
'''Parser for chemical formulas'''

import re
from collections import Counter
import functools
import operator

from .expression.affine import Expression

class FormulaElement(object):
    def __add__(self, other):
        '''Add formula elements creating subformulas'''
        if isinstance(other, FormulaElement):
            if self == other:
                return Formula({ self: 2 })
            return Formula({ self: 1, other: 1 })
        return NotImplemented

    def __radd__(self, other):
        return self + other

    def __or__(self, other):
        '''Merge formula elements into one formula'''
        return Formula({ self: 1 }) | other

    def __ror__(self, other):
        return self | other

    def __mul__(self, other):
        '''Multiply formula element by other'''
        return Formula({ self: other })

    def __rmul__(self, other):
        return self * other

    def repeat(self, count):
        '''Repeat formula element by creating a subformula'''
        return Formula({ self: count })

    def variables(self):
        '''Iterator over variables in formula element'''
        return iter([])

    def substitute(self, **kwargs):
        '''Return formula element with substitutions performed'''
        return self

@functools.total_ordering
class Atom(FormulaElement):
    def __init__(self, symbol):
        self._symbol = symbol

    @property
    def symbol(self):
        return self._symbol

    def __hash__(self):
        return hash('Atom') ^ hash(self._symbol)

    def __eq__(self, other):
        return isinstance(other, Atom) and self._symbol == other._symbol

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        if isinstance(other, Atom):
            return self._symbol < other._symbol
        return NotImplemented

    def __str__(self):
        return self._symbol

    def __repr__(self):
        return 'Atom({})'.format(repr(self._symbol))

class Radical(FormulaElement):
    def __init__(self, symbol):
        self._symbol = symbol

    @property
    def symbol(self):
        return self._symbol

    def __hash__(self):
        return hash('Radical') ^ hash(self._symbol)

    def __eq__(self, other):
        return isinstance(other, Radical) and self._symbol == other._symbol

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        return self._symbol

    def __repr__(self):
        return 'Radical({})'.format(repr(self._symbol))

class Formula(FormulaElement):
    def __init__(self, values={}):
        self._values = {}
        self._variables = set()

        for element, value in values.iteritems():
            if not isinstance(element, FormulaElement):
                raise ValueError('Not a formula element: {}'.format(repr(element)))
            if element != Formula() and value != 0:
                self._values[element] = value

            if callable(getattr(value, 'variables', None)):
                for var in value.variables():
                    self._variables.add(var)
            for var in element.variables():
                self._variables.add(var)

    def substitute(self, **kwargs):
        result = self.__class__()
        for element, value in self._values.iteritems():
            if callable(getattr(value, 'substitute', None)):
                value = value.substitute(**kwargs)
                if isinstance(value, int) and value <= 0:
                    raise ValueError('Expression evaluated to non-positive number')
            result += value * element.substitute(**kwargs)
        return result

    def flattened(self):
        '''Return formula where subformulas have been extracted'''
        stack = [(self, 1)]
        result = {}
        while len(stack) > 0:
            var, value = stack.pop()
            if isinstance(var, Formula):
                for sub_var, sub_value in var._values.iteritems():
                    stack.append((sub_var, value*sub_value))
            else:
                if var in result:
                    result[var] += value
                else:
                    result[var] = value
        return Formula(result)

    def variables(self):
        return iter(self._variables)

    def is_variable(self):
        return len(self._variables) > 0

    def __str__(self):
        '''Return formula represented using Hill notation system'''
        def hill_sorted_elements(values):
            def element_sort_key(pair):
                element, value = pair
                if isinstance(element, Atom):
                    return 0, element.symbol
                elif isinstance(element, Radical):
                    return 1, element.symbol
                else:
                    return 2, None

            if Atom('C') in values:
                yield Atom('C'), values[Atom('C')]
                if Atom('H') in values:
                    yield Atom('H'), values[Atom('H')]
                for element, value in sorted(values.iteritems(), key=element_sort_key):
                    if element not in (Atom('C'), Atom('H')):
                        yield element, value
            else:
                for element, value in sorted(values.iteritems(), key=element_sort_key):
                    yield element, value

        s = ''
        for element, value in hill_sorted_elements(self._values):
            def grouped(element, value):
                return '({}){}'.format(element, value if value != 1 else '')
            def nongrouped(element, value):
                return '{}{}'.format(element, value if value != 1 else '')
            if isinstance(element, Radical):
                s += nongrouped(element, value) if len(element.symbol) == 1 else grouped(element, value)
            elif isinstance(element, Atom):
                s += nongrouped(element, value)
            else:
                s += grouped(element, value)
        return s

    def __repr__(self):
        return 'Formula({})'.format(repr(self._values))

    def __or__(self, other):
        '''Merge formulas into one formula'''
        if isinstance(other, Formula):
            values = Counter(self._values)
            values.update(other._values)
            return Formula(values)
        elif isinstance(other, FormulaElement):
            return self | Formula({ other: 1 })
        return NotImplemented

    def __mul__(self, other):
        '''Multiply formula element by other'''
        values = { key: value*other for key, value in self._values.iteritems() }
        return Formula(values)

    def __hash__(self):
        h = hash('Formula')
        for element, value in self._values.iteritems():
            h ^= hash(element) ^ hash(value)
        return h

    def __eq__(self, other):
        return isinstance(other, Formula) and self._values == other._values

    def __ne__(self, other):
        return not self == other

    @classmethod
    def parse(cls, s):
        '''Parse a formula string (e.g. C6H10O2)'''

        scanner = re.compile(r'''
            (\s+) |         # whitespace
            (\(|\)) |       # group
            ([A-Z][a-z]*) | # element
            (\d+) |         # number
            ([a-z]) |       # variable
            (\Z) |          # end
            (.)             # error
        ''', re.DOTALL | re.VERBOSE)

        def transform_subformula(form):
            '''Extract radical if subformula is a singleton with a radical'''
            if isinstance(form, Formula) and len(form._values) == 1:
                # A radical in a singleton subformula is interpreted as a
                # numbered radical.
                element, value = next(form._values.iteritems())
                if isinstance(element, Radical):
                    return Radical('{}{}'.format(element.symbol, value))
            return form

        stack = []
        formula = Formula()
        expect_count = False

        def close(formula, count=1):
            if len(stack) == 0:
                raise ValueError('Unbalanced parenthesis group in formula')
            subformula = transform_subformula(formula)
            formula = stack.pop()
            if subformula not in formula._values:
                formula._values[subformula] = 0
            formula._values[subformula] += count
            return formula

        for match in re.finditer(scanner, s):
            whitespace, group, element, number, variable, end, error = match.groups()

            if error is not None:
                raise ValueError('Invalid token in formula string: {}'.format(match.group(0)))
            elif whitespace is not None:
                continue
            elif group is not None and group == '(':
                if expect_count:
                    formula = close(formula)
                stack.append(formula)
                formula = Formula()
                expect_count = False
            elif group is not None and group == ')':
                if expect_count:
                    formula = close(formula)
                expect_count = True
            elif element is not None:
                if expect_count:
                    formula = close(formula)
                stack.append(formula)
                if element in 'RX':
                    formula = Radical(element)
                else:
                    formula = Atom(element)
                expect_count = True
            elif number is not None and expect_count:
                formula = close(formula, int(number))
                expect_count = False
            elif variable is not None and expect_count:
                formula = close(formula, Expression(variable))
                expect_count = False
            elif end is not None:
                if expect_count:
                    formula = close(formula)
            else:
                raise ValueError('Invalid token in formula string: {}'.format(match.group(0)))

        if len(stack) > 0:
            raise ValueError('Unbalanced parenthesis group in formula')

        return formula

    @classmethod
    def balance(cls, lhs, rhs):
        '''Return formulas that need to be added to balance given formulas

        Given complete formulas for right side and left side of a reaction,
        calculate formulas for the missing compounds on both sides. Return
        as a left, right tuple. Formulas can be flattened before balancing
        to diregard grouping structure.'''

        def missing(formula, other):
            for element, value in formula._values.iteritems():
                if element not in other._values:
                    yield value*element
                elif value - other._values[element] > 0:
                    delta = value - other._values[element]
                    yield delta*element

        return (reduce(operator.or_, missing(rhs, lhs), Formula()),
                reduce(operator.or_, missing(lhs, rhs), Formula()))

if __name__ == '__main__':
    import doctest
    doctest.testmod()
