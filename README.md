# milp ranker
Derives a ranking from pairwise comparisons (less than or equal, equal, greater than or equal) by solving a mixed-integer linear program (MILP). If the pairwise comparisons contradict themselves (not transitive), the MILP first finds the least changes to find a possible set of transitive comparisons. If not all pairwise comparisons are given, it will return a ranking based on the topological ordering.

Comparisons (i, j) are represented by a number between 0 and 1, where 0 indicates that i <= j, 0.5 that i == j, and 1  that i >= j (and any value in between). `milp_ranker` does not model equal comparisons (`milp_ranker_equal` does model them but it is slower). In `milp_ranker` the cut off between less than and greater than is 0.5 (for `milp_ranker_equal` the 'width' of the equal case can be defined). Setting a comparison to the value of a cut off point, gives the MILP the option to choose the comparison it prefers (it does not model <, >, = it models <=, >=, =).

MILPs are solved with [gurobi](http://www.gurobi.com/).

## Example
Assume we have three nodes we want to compare (A, B, C) and we have the following comparisons: A <= C and B >= C. To use the `milp_ranker`, the node names have to be integers (A=0, B=1, C=2).

```python
from milp_ranker import find_ranking

comparisons = {(0, 2): 0, (1, 2): 1}
ranking, cost = find_ranking(comparisons)
# ranking: [0.0, 2.0, 1.0] -> A has rank 0, B has rank 2, and C has rank 1
# cost: 0.0 -> No changes were made to the comparisons
``` 

If we assume that A <= C and that B == C, we can model the equal case with
```python
from milp_ranker_equal import find_ranking

comparisons = {(0, 2): 0, (1, 2): 0.5}
ranking, cost = find_ranking(comparisons)
# ranking: [0.0, 1.0, 1.0] -> A has rank 0 and B and C share rank 1
# cost: 0.0 -> No changes were made to the comparisons
```

If you would use `milp_ranker` with the equal case, it will make a (random) decision about the ordering:
```python
from milp_ranker import find_ranking

comparisons = {(0, 2): 0, (1, 2): 0.5}
ranking, cost = find_ranking(comparisons)
# ranking: [0.0, 0.0, 1.0] -> A and B have rank 0 and C has rank 1
# cost: 0.0 -> No changes were made to the comparisons
```
