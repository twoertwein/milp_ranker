import numpy as np
from gurobipy import GRB, LinExpr, Model


def find_ranking(comparisons, equal_width=0.2, max_rank=-1):
    """
    Find minimal changes to a set of comparisons so that they are consistent
    (transitive), it return a topological ranking.

    comparisons     A dictionary with tuple keys in the form of (i, j), values
                    are scalars indicating the probability of i > j. It is
                    assumed that comparisons are symmetric. Use 0 for i < j,
                    0.5 for i == j, and 1 for i > j (and any value in between).
    equal_width     0..0.5-equal_width/2 is considered '<=' and 0.5..0.5+equal_width/2
                    is considered '>='. In between it is considered to be '=='.
    max_rank        Maximal rank, a low value forces the model to choose more
                    equal cases.

    Returns:
        A tuple of size two:
        0) Ranking derived from topological sort (list of ranks in order of
           nodes);
        1) Sum of absolute changes to the comparisons.
    """
    # remove unnecessary variables
    comparisons = {(i, j) if i < j else (j, i): value if i < j else 1 - value
                   for (i, j), value in comparisons.items()}
    nodes = np.unique(
        [i for ij in comparisons.keys() for i in ij])

    # define variables
    model = Model('comparison')
    values = np.fromiter(comparisons.values(), dtype=float)
    assert values.max() <= 1 and values.min() >= 0
    # variables to encode the error of comparisons
    E_ij = model.addVars(comparisons.keys(), name='e_ij', vtype=GRB.CONTINUOUS,
                         ub=1.0-values, lb=-values)
    # variables to encode hard choice of >=, <=, ==
    Ge_ij = model.addVars(comparisons.keys(), name='ge_ij', vtype=GRB.BINARY)
    Le_ij = model.addVars(comparisons.keys(), name='le_ij', vtype=GRB.BINARY)
    Eq_ij = model.addVars(comparisons.keys(), name='eq_ij', vtype=GRB.BINARY)
    # variables to help with transitivity in non-fully connected graphs
    if max_rank < 1:
        max_rank = len(nodes)
    R_i = model.addVars(nodes, name='r_i', vtype=GRB.CONTINUOUS, lb=0,
                        ub=max_rank)
    # variables to emulate abs
    T_ij_pos = {}
    T_ij_neg = {}
    index = (values != 1) & (values != 0)
    T_ij_pos = model.addVars(
        (ij for ij, value in comparisons.items() if value not in [0.0, 1.0]),
        vtype=GRB.CONTINUOUS, name='T_ij_pos', lb=0, ub=1-values[index])
    T_ij_neg = model.addVars(
        (ij for ij, value in comparisons.items() if value not in [0.0, 1.0]),
        vtype=GRB.CONTINUOUS, name='T_ij_neg', lb=0, ub=values[index])
    model.update()

    # emulate abs for non-binary comparisons: E_ij = T_ij_pos - T_ij_neg
    model.addConstrs(
        (E_ij[ij] == T_ij_pos[ij] - T_ij_neg[ij] for ij in T_ij_pos),
        'E_ij = T_ij_pos - T_ij_neg')

    # hard decision of >=, <=, and ==
    lower_bound = 0.5 - equal_width / 2.0
    upper_bound = 0.5 + equal_width / 2.0
    # <=
    model.addConstrs(
        (E_ij[ij] + comparisons[ij] - upper_bound <= ge_ij
         for ij, ge_ij in Ge_ij.items()), 'ge_ij_lower_bound')
    model.addConstrs(
        (E_ij[ij] + comparisons[ij] - upper_bound >= -1 + ge_ij
         for ij, ge_ij in Ge_ij.items()), 'ge_ij_upper_bound')
    # >=
    model.addConstrs(
        (E_ij[ij] + comparisons[ij] - lower_bound >= -le_ij
         for ij, le_ij in Le_ij.items()), 'le_ij_lower_bound')
    model.addConstrs(
        (E_ij[ij] + comparisons[ij] - lower_bound <= 1 - le_ij
         for ij, le_ij in Le_ij.items()), 'le_ij_upper_bound')
    # ==
    model.addConstrs(
        (le + eq + ge == 1 for le, eq, ge in zip(
            Le_ij.values(), Eq_ij.values(), Ge_ij.values())), 'eq_ij')

    # transitivity
    for (i, j), eq_a in Eq_ij.items():
        le_a = Le_ij[i, j]
        ge_a = Ge_ij[i, j]
        for k in nodes:
            j_, k_ = j, k
            if j > k:
                j_, k_ = k, j
            eq_b = Eq_ij.get((j_, k_), None)
            if eq_b is None:
                continue
            else:
                le_b = Le_ij[j_, k_]
                ge_b = Ge_ij[j_, k_]
            if j_ != j:
                le_b, ge_b = ge_b, le_b

            i_, k_ = i, k
            if i > k:
                i_, k_ = k, i
            eq_c = Eq_ij.get((i_, k_), None)
            if eq_c is None:
                continue
            else:
                le_c = Le_ij[i_, k_]
                ge_c = Ge_ij[i_, k_]
            if i_ != i:
                le_c, ge_c = ge_c, le_c

            # a <= b and b <= c -> a <= c
            model.addLConstr(
                ge_a + ge_b, GRB.LESS_EQUAL, 1 + ge_c,
                f'transitivity_ge_{i},{j},{k}')
            # a >= b and b >= c -> a >= c
            model.addLConstr(
                le_a + le_b, GRB.LESS_EQUAL, 1 + le_c,
                f'transitivity_le_{i},{j},{k}')
            # a <= b and b == c -> a <= c
            model.addLConstr(
                le_a + eq_b, GRB.LESS_EQUAL, 1 + le_c,
                f'transitivity_leeq_{i},{j},{k}')
            # a == b and b <= c -> a <= c
            model.addLConstr(
                eq_a + le_b, GRB.LESS_EQUAL, 1 + le_c,
                f'transitivity_eqle_{i},{j},{k}')
            # a >= b and b == c --> a >= c
            model.addLConstr(
                ge_a + eq_b, GRB.LESS_EQUAL, 1 + ge_c,
                f'transitivity_geeq_{i},{j},{k}')
            # a == b and b >= c --> a >= c
            model.addLConstr(
                eq_a + ge_b, GRB.LESS_EQUAL, 1 + ge_c,
                f'transitivity_eqge_{i},{j},{k}')
            # a == b and b == c --> a == c
            model.addLConstr(
                eq_a + eq_b, GRB.LESS_EQUAL, 1 + eq_c,
                f'transitivity_eq_{i},{j},{k}')

    # transitivity helper (for not-fully connected graphs)
    # also provides a latent rank
    big_m = max_rank
    model.addConstrs(
        ((1 - ge_ij) * big_m + R_i[i] >= R_i[j] + 1 for (i, j), ge_ij in Ge_ij.items()),
        'rank_transitivity_larger')
    model.addConstrs(
        ((1 - le_ij) * big_m + R_i[j] >= R_i[i] + 1 for (i, j), le_ij in Le_ij.items()),
        'rank_transitivity_smaller')
    model.addConstrs(
        ((1 - eq_ij) * big_m + R_i[j] >= R_i[i] for (i, j), eq_ij in Eq_ij.items()),
        'rank_transitivity_equal1')
    model.addConstrs(
        ((1 - eq_ij) * big_m + R_i[i] >= R_i[j] for (i, j), eq_ij in Eq_ij.items()),
        'rank_transitivity_equal2')

    # objective function
    objective = LinExpr()
    for ij, value in comparisons.items():
        if value == 1.0:
            objective += -E_ij[ij]
        elif value == 0.0:
            objective += E_ij[ij]
        else:
            objective += T_ij_pos[ij] + T_ij_neg[ij]
    model.setObjective(objective, GRB.MINIMIZE)

    # solve
    model.optimize()

    # verify abs emulation: one T_ij has to be 0
    for ij, value in T_ij_pos.items():
        assert value.X == 0 or T_ij_neg[ij] == 0, \
            f'T_{ij} pos {value.X} neg {T_ij_neg[ij]}'

    # find minimal Rs
    model_ = Model('comparison')
    R_i = model_.addVars(nodes, name='r_i', vtype=GRB.CONTINUOUS, lb=0,
                         ub=len(nodes))
    for ((i, j), ge_ij), le_ij in zip(Ge_ij.items(), Le_ij.values()):
        if ge_ij.x == 1:
            model_.addConstr(R_i[i] >= R_i[j] + 1)
        elif le_ij.x == 1:
            model_.addConstr(R_i[j] >= R_i[i] + 1)
        else:
            model_.addConstr(R_i[j] == R_i[i])
    model_.setObjective(R_i.sum(), GRB.MINIMIZE)
    model_.optimize()

    return [model_.getVarByName(f'r_i[{i}]').X for i in range(len(nodes))], \
        model.getObjective().getValue()
