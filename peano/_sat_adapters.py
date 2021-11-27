"""
This module provides one class, CurveSATAdapter,
which allows to encode fuzzy curves by boolean formulas
and apply SAT-solvers to find a curve satisfying clauses.
"""
import itertools
from collections import namedtuple

from pysat.solvers import Glucose3

from .curves import Curve


class CurveSATAdapter:
    """
    Class-interface between fuzzy curves and SAT-solvers.

    Fuzzy curve is encoded by boolean variables:
    For each fraction (pnum, cnum) and each allowed spec S
    there is a boolean variable x s.t. x=True <=> spec at (pnum,cnum) is S
    The crucial property of a fuzzy curve is that specs for
    different position are independent of each other (so we do not need
    additional clauses), this is satisfied for PathFuzzyCurve class.

    Restrictions on curve specs (clauses) are written and CNF formulas.
    SAT-Solver can find a model, i.e. curve satisfying clauses.

    We maintain that CNF formulas, each clause (disjunction) is represented by dict {var: True|False}
    var is a hashable object, there are several types of variables:
    - sp_var: a variable for each variant of (pnum, cnum, spec)
    - junc_var: a variable for each junction
    - curve_var: variables for some fuzzy curves (internal)
    """

    def __init__(self, curve=None):
        """
        Create SAT adapter for given fuzzy curve.

        Inits basic vars and clauses.
        """
        self._int_clauses = []
        self._curve_vars = set()  # just cache (?)
        self._var_no = {}
        self.solver = None  # set in solve()

        if curve is not None:
            self._load_curve(curve)

    def _load_curve(self, curve):
        self.curve = curve  # used in get_model

        # possible specs
        for pnum in range(curve.pcount):
            for cnum in range(curve.genus):
                self._make_only([self._get_sp_var(pnum, cnum, sp) for sp in curve.gen_allowed_specs(pnum, cnum)])

        # auto-junction - in each curve
        for junc in curve.gen_auto_junctions():
            self._append_clause({self._get_junc_var(junc): True})

        # regular junctions
        for junc, curves in curve.get_junction_templates().items():
            self._make_junc_var(junc, curves)

    def copy(self):
        assert self.solver is None
        new_adapter = type(self)()
        new_adapter._int_clauses = self._int_clauses[:]
        new_adapter._curve_vars = self._curve_vars.copy()
        new_adapter._var_no = self._var_no.copy()
        new_adapter.curve = self.curve
        return new_adapter

    @staticmethod
    def _get_sp_var(pnum, cnum, sp):
        # sp_var: var=True <=> curve has given spec at given pnum, cnum.
        return 'spec', pnum, cnum, sp

    @staticmethod
    def _get_junc_var(junc):
        # See _make_junc_var."""
        return 'junc', junc

    def _make_only(self, only_vars):
        # Add clauses that one and only one var is True.
        # (var1 or var2 or ... or varN) and (!var_i or !var_j)
        # N.B. try to use add_atmost
        self._append_clause({var: True for var in only_vars})
        for var_pair in itertools.combinations(only_vars, 2):
            self._append_clause({var_pair[0]: False, var_pair[1]: False})

    def _make_curve_var(self, curve):
        # Create variable Z, such that Z=True <=> maintained curve is consistent with given one.
        curve_info = tuple(curve.gen_defined_specs())
        Z = ('curve', curve_info)
        if Z in self._curve_vars:
            return Z  # already initiated

        # Z <-> curve  <=>  (Z->curve) and (curve->Z)
        # curve <=> (sp1 and sp2 ... and spk)

        # Z->curve  <=>  !Z or curve  <=>  (!Z or sp1) and (!Z or sp2) and ... (!Z or spk)
        for pnum, cnum, sp in curve_info:
            self._append_clause({Z: False, self._get_sp_var(pnum, cnum, sp): True})

        # curve->Z  <=>  !curve or Z  <=>  !sp1 or !sp2 or ... or !spk or Z
        clause_rev = {self._get_sp_var(pnum, cnum, sp): False for pnum, cnum, sp in curve_info}
        clause_rev[Z] = True
        self._append_clause(clause_rev)

        self._curve_vars.add(Z)
        return Z
        
    def _make_junc_var(self, junc, curves):
        # Create variable J, such that J=True <=> curve has junction J
        # curves is the list of minimal fuzzy curves with this junc (see get_junction_templates)
        J = self._get_junc_var(junc)

        curve_vars = [self._make_curve_var(curve) for curve in curves]
        # J <-> (c1 or c2 or .. or ck)

        # J->(c1 or c2 or .. or ck)  <=>  !J or c1 or c2 .. or ck
        clause = {cv: True for cv in curve_vars}
        clause[J] = False
        self._append_clause(clause)

        # (c1 or c2 .. or ck)->J  <=>  !(c1 or .. or ck) or J  <=>  (!c1 and ... !ck) or J
        #                         <=>  (!c1 or J) and ... (!ck or J)
        for cv in curve_vars:
            self._append_clause({cv: False, J: True})

        return J

    def add_spec_clause(self, pnum, cnum, spec):
        """Demand that maintained fuzzy curve has given spec."""
        self._append_clause({self._get_sp_var(pnum, cnum, spec): True})

    def add_forbid_clause(self, junc, curve):
        """Forbid that maintained fuzzy curve has given junction and is consistent with given sub-curve."""
        # !(J and bm1 and bm2 .. and bmk) = !J or !bm1 or !bm1 ..
        junc_var = self._get_junc_var(junc)
        clause = {self._get_sp_var(pnum, cnum, sp): False for pnum, cnum, sp in curve.gen_defined_specs()}
        clause[junc_var] = False
        self._append_clause(clause)

    def _append_clause(self, clause):
        # Add clause to CNF, maintain int representation."""
        int_clause = []
        for var, val in clause.items():
            if var not in self._var_no:
                max_var_no = 1 + len(self._var_no)
                self._var_no[var] = max_var_no
            var_no = self._var_no[var]
            token = var_no if val else -var_no
            int_clause.append(token)
        self._int_clauses.append(tuple(int_clause))

    SATProblemSize = namedtuple('SATProblemSize', ['literals', 'clauses', 'variables'])

    def get_problem_size(self):
        return self.SATProblemSize(
            literals=sum(len(clause) for clause in self._int_clauses),
            clauses=len(self._int_clauses),
            variables=len(self._var_no),
        )

    def solve(self):
        """
        Call SAT-solver to determine if there is a curve satisfying clauses.

        We create Solver object from scratch every time.

        Returns:
            True if there is a model, False otherwise
        """
        self.solver = Glucose3()
        self.solver.append_formula(self._int_clauses)
        return self.solver.solve()

    def get_model(self):
        """
        Get a model, i.e., values for variables satisfying clauses.

        Only active variables are used in model, for other vars
        any value may be set.

        A noteworthy restriction of used SAT-solver is that it
        returns only one model, not the list of all possible models.

        Returns:
            dict {var: value}
        """
        int_model = self.solver.get_model()
        model = {}
        no_var = {var_no: var for var, var_no in self._var_no.items()}
        for int_tok in int_model:
            var = no_var[abs(int_tok)]
            model[var] = (int_tok > 0)
        return model

    def get_curve_from_model(self, model):
        """
        Get a full-defined curve from a model.

        Returns:
            Curve instance

        Raises:
            ValueError: if can't find a curve (bad model)
        """
        # since curve is encoded by spec_vars and other vars
        # are expressed in terms of them (see, e.g., _make_junc_var),
        # we should check only spec_vars
        curve = self.curve
        specs_allowed_by_model = []
        G = curve.genus
        for pnum in range(curve.pcount):
            for cnum in range(G):
                good_spec = None
                for sp in curve.gen_allowed_specs(pnum, cnum):
                    sp_var = self._get_sp_var(pnum, cnum, sp)
                    if (sp_var not in model) or model[sp_var]:
                        good_spec = sp
                        break
                if good_spec is None:
                    raise ValueError("Bad model, can't find curve!")
                specs_allowed_by_model.append(sp)

        patterns = []
        for pnum, pattern in enumerate(curve.patterns):
            specs = specs_allowed_by_model[pnum * G : (pnum + 1) * G]
            patterns.append((pattern.proto, specs))
        return Curve(dim=curve.dim, div=curve.div, patterns=patterns)
