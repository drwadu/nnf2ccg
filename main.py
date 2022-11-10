from clingo import Control
from typing import Dict, List
from functools import reduce
from operator import add, mul
from math import log10


AND: int = 1
OR: int = 0


def transpile(
    nnf_path: str, lp_path: str, cnf_path: str, verbose: bool = False
) -> None:
    name = nnf_path.split(".")[0]

    lp = open(lp_path, "r").read()
    ctl = Control("0")
    ctl.add("base", [], lp)
    ctl.ground([("base", [])])

    cnf_mappings: Dict[str, int] = {
        atom_str: int(atom_int)
        for (atom_int, atom_str) in map(
            lambda l: l.rstrip("\n").split(" ")[1:],
            filter(lambda l: l.startswith("c "), open(cnf_path, "r").readlines()),
        )
    }
    atom_mappings, supported_atoms = dict(), []
    for atom in ctl.symbolic_atoms:  # one traversal vs. comprehension
        s = str(atom.symbol)
        atom_mappings[cnf_mappings[s]] = atom.literal
        supported_atoms.append(s)
    falsified = list(
        map(
            lambda a: cnf_mappings[a],
            filter(lambda a: a not in supported_atoms, cnf_mappings.keys()),
        )
    )

    nnf = map(lambda l: l.rstrip("\n"), open(nnf_path, "r").readlines().__iter__())
    stats = next(nnf, None)
    node_count = int(stats.split(" ")[1])

    nodes, n_nodes, n_edges = [], 0, 0

    atom_popped_ids, popped_ids = set(), set()
    new_vars_count, n_popped = 0, 0
    node_id_diffs: Dict[int, int] = dict()

    line, i = next(nnf, None), 0
    while line:
        spec = line.split(" ").__iter__()
        node_type = next(spec, None)
        node = []

        if node_type == "L":
            lit = int(next(spec, None))
            atom = abs(lit)

            if not atom_mappings.get(atom, None) or atom in falsified:
                nodes.append([])
                atom_popped_ids.add(i)
                n_popped += 1
            else:
                n_nodes += 1
                if lit > 0:
                    new_vars_count += 1
                nodes.append([1, lit])
                node_id_diffs[i] = i - n_popped
        else:
            is_or_node = node_type == "O"

            if is_or_node:
                next(spec, None)
            next(spec, None)

            children = list(
                map(
                    lambda t: (t[0], nodes[t[0]], t[1]),
                    filter(
                        lambda t: t[1] < n_nodes,
                        map(
                            lambda c: (c, node_id_diffs.get(c, c)),
                            filter(lambda c: c not in atom_popped_ids, map(int, spec)),
                        ),
                    ),
                )
            )
            n_children = len(children)

            if n_children == 0:
                popped_ids.add(i)
                node_id_diffs[i] = node_count
                n_popped += 1
            elif n_children == 1:
                is_root = i == node_count - 1
                if not is_root:
                    _, nnf_child_node, ccg_child_id = children[0]
                    node = nnf_child_node
                    node_id_diffs[i] = ccg_child_id
                n_popped += 1
                popped_ids.add(i)
            else:
                n_nodes += 1
                n_edges += n_children

                if not is_or_node:
                    val = reduce(mul, map(lambda t: next(iter(t[1]), 1), children))
                    node.append(val)
                    node.extend(map(lambda t: t[2], children))
                    node.append(AND)
                else:
                    val = reduce(add, map(lambda t: next(iter(t[1]), 0), children))
                    node.append(val)
                    node.extend(map(lambda t: t[2], children))
                    node.append(OR)

                node_id_diffs[i] = i - n_popped
            nodes.append(node)

        line = next(nnf, None)
        i += 1

    ccg = list(
        map(lambda t: t[1], 
        filter(
            lambda t: t[0] not in atom_popped_ids and t[0] not in popped_ids,
            enumerate(nodes),
        ))
    )

    node_count_ccg = len(ccg)
    assert node_count_ccg == node_count - n_popped
    assert node_count_ccg == n_nodes

    root = ccg[node_count_ccg - 1]
    count = root[0]

    stats = f"ccg {node_count_ccg} {n_edges} {new_vars_count} {log10(count)}"

    print(stats)
    #map(lambda t: print(f"c {t[0]} {t[1]}"), cnf_mappings)
    for a,b in cnf_mappings.items():
        print(f"c {b} {a}")
    for node in ccg:
        node_size = len(node)
        if node_size == 2:
            for s in reversed(node):
                print(f"{s}", end=" ")
            print()
        else:
            children_count = node_size - 2
            last_idx = node_size - 1
            node_type = "*" if node[last_idx] == AND else "+"
            print(node_type, end=" ")
            print(children_count, end=" ")
            for s in node[1:last_idx]:
                print(f"{s}", end=" ")
            print()



if __name__ == "__main__":
    import sys

    nnf_path = sys.argv[1]
    lp_path = sys.argv[2]
    cnf_path = sys.argv[3]

    transpile(nnf_path, lp_path, cnf_path)
