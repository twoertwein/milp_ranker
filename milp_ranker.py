import numpy as np
from gurobipy import GRB, LinExpr, Model


def find_ranking(comparisons):
    """
    Find minimal changes to a set of comparisons so that they are consistent
    (transitive), it return a topological ranking.

    comparisons     A dictionary with tuple keys in the form of (i, j), values
                    are scalars indicating the probability of i >= j. It is
                    assumed that comparisons are symmetric. The only way this
                    function it can handle '=' is by setting the probability to
                    0.5. In this case the MILP can treat it as > or < (but not
                    as =).

    Returns:
        A tuple of size two:
        0) Ranking derived from topological sort (list of ranks in order of nodes);
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
    # variable to encode hard choice of >= and <=
    Z_ij = model.addVars(comparisons.keys(), name='z_ij', vtype=GRB.BINARY)
    # variables to help with transitivity in non-fully connected graphs
    R_i = model.addVars(nodes, name='r_i', vtype=GRB.CONTINUOUS, lb=0,
                        ub=len(nodes))
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

    # hard decision of >= and <=: z_ij == 1 <-> i > j
    model.addConstrs(
        (E_ij[ij] + comparisons[ij] - 0.5 >= - 1 + z_ij
         for ij, z_ij in Z_ij.items()), 'z_ij_upper_bound')
    model.addConstrs(
        (E_ij[ij] + comparisons[ij] - 0.5 <= z_ij
         for ij, z_ij in Z_ij.items()), 'z_ij_lower_bound')

    # transitivity
    for (i, j), a in Z_ij.items():
        for k in nodes:
            j_, k_ = j, k
            if j > k:
                j_, k_ = k, j
            b = Z_ij.get((j_, k_), None)
            if b is None:
                continue
            elif j_ != j:
                b = 1 - b

            i_, k_ = i, k
            if i > k:
                i_, k_ = k, i
            c = Z_ij.get((i_, k_), None)
            if c is None:
                continue
            elif i_ != i:
                c = 1 - c
            # a <= b and b <= c -> a <= c
            model.addLConstr(
                a + b, GRB.LESS_EQUAL, 1 + c,
                f'transitivity_ge_{i},{j},{k}')
            # a >= b and b >= c -> a >= c
            model.addLConstr(
                a + b, GRB.GREATER_EQUAL, c,
                f'transitivity_le_{i},{j},{k}')

    # transitivity helper (for not-fully connected graphs)
    # also provides a latent rank
    big_m = len(nodes)
    model.addConstrs(
        ((1 - z_ij) * big_m + R_i[i] >= R_i[j] + 1 for (i, j), z_ij in Z_ij.items()),
        'rank_transitivity_larger')
    model.addConstrs(
        (z_ij * big_m + R_i[j] >= R_i[i] + 1 for (i, j), z_ij in Z_ij.items()),
        'rank_transitivity_smaller')

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
    for (i, j), z_ij in Z_ij.items():
        if z_ij.x == 1:
            model_.addConstr(R_i[i] >= R_i[j] + 1)
        else:
            model_.addConstr(R_i[j] >= R_i[i] + 1)
    model_.setObjective(R_i.sum(), GRB.MINIMIZE)
    model_.optimize()

    return [model_.getVarByName(f'r_i[{i}]').X for i in range(len(nodes))], \
        model.getObjective().getValue()
