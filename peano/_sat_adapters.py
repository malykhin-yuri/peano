import itertools
from collections import namedtuple

from pysat.solvers import Glucose3

from .curves import Curve


class CurveSATAdapter:
    # Class-interface between fuzzy curves and SAT-solvers.

    # Fuzzy curve is encoded by boolean variables:
    # For each fraction (pnum, cnum) and each allowed spec S
    # there is a boolean variable x s.t. x=True <=> spec at (pnum,cnum) is S
    # It is required that that specs for # different positions
    # are independent of each other, so we do not need
    # additional clauses. This is satisfied for the PathFuzzyCurve class.

    # Restrictions on curve specs (clauses) are written as CNF formulas.
    # SAT-Solver can find a model, i.e. curve satisfying clauses.

    # We maintain arising CNF formulas, each clause (disjunction) uses variables:
    # - sp_var: base variables, for each variant of (pnum, cnum, spec)
    # - junc_var: additional variable for to represent presence of a junction
    # - curve_var: additional variables for some fuzzy curves (internal)

    def __init__(self, curve=None, adapter=None):
        self._int_clauses = []
        self._curve_vars = set()  # just cache (?)
        self._var_no = {}
        self._solver = None  # set in solve()

        if curve is not None:
            self._load_curve(curve)
        elif adapter is not None:
            self._copy_from(adapter)
        else:
            raise ValueError("Init with curve or adapter!")

    def _load_curve(self, curve):
        self._curve = curve  # used in get_model

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

    def _copy_from(self, other):
        assert other._solver is None
        self._int_clauses = other._int_clauses[:]
        self._curve_vars = other._curve_vars.copy()
        self._var_no = other._var_no.copy()
        self._curve = other._curve

    @staticmethod
    def _get_sp_var(pnum, cnum, sp):
        # sp_var: var=True <=> curve has given spec at given pnum, cnum.
        return 'spec', pnum, cnum, sp

    @staticmethod
    def _get_junc_var(junc):
        # See _make_junc_var
        return 'junc', junc

    @staticmethod
    def _get_curve_var(curve):
        # See _make_curve_var
        return 'curve', tuple(curve.gen_defined_specs())

    def _make_only(self, only_vars):
        # Add clauses that exactly one var is True.
        # (var1 or var2 or ... or varN) and (!var_i or !var_j)
        # N.B. try to use add_atmost
        self._append_clause({var: True for var in only_vars})
        for var1, var2 in itertools.combinations(only_vars, 2):
            self._append_clause({var1: False, var2: False})

    def _make_curve_var(self, curve):
        # Create variable Z, such that Z=True <=> maintained curve is consistent with given one.
        Z = self._get_curve_var(curve)
        if Z in self._curve_vars:
            return Z  # already initiated

        curve_info = list(curve.gen_defined_specs())
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

    def add_forbid_clause(self, junc, curve):
        # Forbid that maintained fuzzy curve has given junction and is consistent with given sub-curve
        # !(J and sp1 and sp2 .. and spk) = !J or !sp1 or !sp1 ..
        junc_var = self._get_junc_var(junc)
        clause = {self._get_sp_var(pnum, cnum, sp): False for pnum, cnum, sp in curve.gen_defined_specs()}
        clause[junc_var] = False
        self._append_clause(clause)

    def _append_clause(self, clause):
        # Add clause to CNF, maintain int representation
        int_clause = []
        for var, value in clause.items():
            if var not in self._var_no:
                max_var_no = 1 + len(self._var_no)
                self._var_no[var] = max_var_no
            var_no = self._var_no[var]
            token = var_no if value else -var_no
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
        # Call SAT-solver to determine if there is a curve satisfying clauses.
        # We create Solver object from scratch every time.
        # Returns:
        #     True if there is a model, False otherwise
        self._solver = Glucose3()
        self._solver.append_formula(self._int_clauses)
        return self._solver.solve()

    def _get_model(self):
        # Get a model, i.e., values for variables satisfying clauses.
        # Only active variables are used in model, for other vars any value may be set.
        # A noteworthy restriction of used SAT-solver is that it
        # returns only one model, not the list of all possible models.
        int_model = self._solver.get_model()
        model = {}
        no_var = {var_no: var for var, var_no in self._var_no.items()}
        for int_tok in int_model:
            var = no_var[abs(int_tok)]
            model[var] = (int_tok > 0)
        return model

    def get_model_curve(self):
        # Get a full-defined curve corresponding to found model
        model = self._get_model()

        # since curve is encoded by spec_vars and other vars
        # are expressed in terms of them (see, e.g., _make_junc_var),
        # we should check only spec_vars
        curve = self._curve
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
