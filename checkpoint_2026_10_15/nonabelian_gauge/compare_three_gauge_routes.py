#!/usr/bin/env python3
"""Reproduce the cross-comparison of the independently authored UV routes.

The two integration cores are standalone and do not import each other. This
separate comparison harness deliberately imports both adjacent core files.
Only the Python standard library is used. No files are written, and explicit
failures remain active under python -O.
"""
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from pathlib import Path
import json


def load_adjacent(name):
    path = Path(__file__).resolve().parent / (name + '.py')
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError('cannot load adjacent comparison module: ' + name)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def independent_name(label):
    if label[0] == 'O':
        return 'O_F' + ''.join(map(str, label[1])) + '_S' + ''.join(map(str, label[2]))
    if label[0] == 'GF':
        return 'SG_F' + str(label[2]) + '_S' + ''.join(map(str, label[1]))
    if label[0] == 'GS':
        if label[3] == 'GthenS':
            blocks = (label[1], (label[2],))
        elif label[3] == 'SthenG':
            blocks = ((label[2],), label[1])
        else:
            raise RuntimeError('unknown scalar-seagull ordering')
        return 'SG_F_S' + '_'.join(''.join(map(str, block)) for block in blocks)
    raise RuntimeError('unknown producer topology label')


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def main():
    producer = load_adjacent('verify_three_gauge')
    auditor = load_adjacent('verify_three_gauge_cubature')
    independent_topologies = {top['name']: top for top in auditor.TOPS}
    require(len(independent_topologies) == 33, 'expected 33 distinct independent topologies')
    require(producer.B == auditor.BASIS, 'the polarization bases differ')
    individual_comparisons = 0
    tensor_comparisons = 0
    topology_mapping = None
    for components in product(range(4), repeat=3):
        eps = tuple(auditor.BASIS[j] for j in components)
        diagrams = producer.raw_diagrams(eps)
        require(len(diagrams) == 33, 'expected 33 producer topologies')
        mapped = [independent_name(label) for label, _ in diagrams]
        require(len(set(mapped)) == 33, 'topology mapping is not one-to-one')
        require(set(mapped) == set(independent_topologies), 'topology mapping is incomplete')
        if topology_mapping is None:
            topology_mapping = mapped
        else:
            require(topology_mapping == mapped, 'topology mapping changed between polarizations')
        for (label, expected), name in zip(diagrams, mapped):
            actual = auditor.graph(independent_topologies[name], eps)
            require(actual == expected, 'individual graph mismatch at ' + str((components, label)))
            individual_comparisons += 1
        actual_words, _ = auditor.integrate(eps)
        actual_tensor = {(word, mask): coefficient
                         for word, part in actual_words.items()
                         for mask, coefficient in part.items()}
        expected_tensor = producer.colored_total(eps)
        require(actual_tensor == expected_tensor, 'free-word tensor mismatch at ' + str(components))
        require(all(isinstance(value, (int, Fraction)) for value in actual_tensor.values()),
                'unexpected non-exact tensor coefficient')
        tensor_comparisons += 1
    require(individual_comparisons == 2112, 'incomplete individual graph comparison')
    require(tensor_comparisons == 64, 'incomplete complete-tensor comparison')
    result = {
        'status': 'PASS',
        'independent_cubature_vs_producer_pairings_individual_graph_comparisons': individual_comparisons,
        'full_free_word_tensor_basis_comparisons': tensor_comparisons,
        'producer': 'verify_three_gauge.py',
        'independent_auditor': 'verify_three_gauge_cubature.py',
        'comparison_harness': 'compare_three_gauge_routes.py',
        'comparison_imports_both_cores': True,
        'core_routes_remain_standalone': True,
        'arithmetic': 'exact standard-library integers and Fractions; cubature uses Gaussian integer pairs',
        'meaning': 'Every one of the 33 graph residues agrees at all 64 polarization basis triples before the separate 6-free-color-word total comparison.'
    }
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
