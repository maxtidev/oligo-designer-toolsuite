import os
from collections import defaultdict
from itertools import combinations
from typing import Literal

import yaml

OUTPUT_PATH = "validation/outputs"

BASE = "baseline"
OTHER = "dp"

WORST_RELEVANT = {"merfish", "scrinshot", "seqfishplus"}
AVERAGE_RELEVANT = {"oligo", "cyclehcr", "hcr"}


def read_oligosets_paths(subpath: str) -> dict[str, str]:
    # returns pipeline_name -> yaml path
    result = {}

    base_path = os.path.join(OUTPUT_PATH, subpath)
    for pipeline_output_dirname in os.listdir(base_path):
        pipeline_name = pipeline_output_dirname.split("_")[1]

        pipeline_output_path = os.path.join(base_path, pipeline_output_dirname)
        for output_filename in os.listdir(pipeline_output_path):
            output_file_path = os.path.join(pipeline_output_path, output_filename)
            if output_filename.endswith("_probes.yml"):
                result[pipeline_name] = output_file_path
                break

    return result


def matching_score_tuple(
    scores_dict: dict[str, float], relevant: Literal["worst", "average"]
) -> tuple[float, float]:
    if relevant == "worst":
        return (scores_dict["set_score_worst"], scores_dict["set_score_sum"])
    else:
        return (scores_dict["set_score_average"], scores_dict["set_score_worst"])


def extract_scores(
    region_oligosets_dict: dict, relevant: Literal["worst", "average"]
) -> dict[str, list[dict[str, float]]]:
    # returns region_id -> [score_type -> score (set 1), ..., score_type -> score (set n)]
    result: defaultdict[str, list[dict[str, float]]] = defaultdict(list)

    for region_id, oligosets_dict in region_oligosets_dict.items():
        for oligoset_dict in oligosets_dict.values():
            scores = oligoset_dict["Oligoset Score"]
            result[region_id].append(scores)

        # normalize order
        # remember that lower scores are better, so this sorts sets best-to-worst
        result[region_id].sort(
            key=lambda scores: matching_score_tuple(scores, relevant),
        )

    return result


def get_jaccard_stats(oligosets_dict: dict) -> tuple[float, float]:

    oligosets: list[set[str]] = []

    for oligoset_dict in oligosets_dict.values():
        oligoset = set()
        for key, oligo_dict in oligoset_dict.items():
            if key == "Oligoset Score":
                continue
            oligoset.add(oligo_dict["oligo_id"])
        oligosets.append(oligoset)

    if len(oligosets) < 2:
        return (0, 0)

    jaccard_indices = []
    for oligoset_1, oligoset_2 in combinations(oligosets, 2):
        jaccard_index = len(oligoset_1 & oligoset_2) / len(oligoset_1 | oligoset_2)
        jaccard_indices.append(jaccard_index)

    return (max(jaccard_indices), sum(jaccard_indices) / len(jaccard_indices))


def compare_oligosets(
    base_oligosets: list[dict[str, float]],
    other_oligosets: list[dict[str, float]],
    relevant: Literal["worst", "average"],
) -> None:
    is_other_better: bool = True

    print("Length:", len(base_oligosets), len(other_oligosets), end=" ")

    if len(other_oligosets) < len(base_oligosets):
        is_other_better = False
        print("❌❌❌")
    else:
        print("✔️")

    for base_scores, other_scores in zip(base_oligosets, other_oligosets):
        base_scores_tuple = matching_score_tuple(base_scores, relevant)
        other_scores_tuple = matching_score_tuple(other_scores, relevant)

        print(base_scores_tuple, other_scores_tuple, end=" ")

        # remember that lower scores are better
        if other_scores_tuple > base_scores_tuple:
            is_other_better = False
            print("❌")
        else:
            print("✔️")

    print("Is other better:", "YES ✔️" if is_other_better else "NO ❌")


def compare_all_oligosets(
    base_yml_path: str, other_yml_path: str, relevant: Literal["worst", "average"]
) -> None:
    # read into dict, sort keys and value dicts, compare
    with open(base_yml_path) as f:
        base_dict: dict = yaml.safe_load(f)
    with open(other_yml_path) as f:
        other_dict: dict = yaml.safe_load(f)

    base_scores = extract_scores(base_dict, relevant)
    other_scores = extract_scores(other_dict, relevant)

    for region_id, base_oligosets in base_scores.items():
        print(f"Comparing {region_id}...")

        base_jaccard_idx_max, base_jaccard_idx_avg = get_jaccard_stats(base_dict[region_id])
        other_jaccard_idx_max, other_jaccard_idx_avg = get_jaccard_stats(other_dict[region_id])

        print(f"Highest similarity: {base_jaccard_idx_max:.2f} {other_jaccard_idx_max:.2f}", end=" ")
        print("↗️" if other_jaccard_idx_max > base_jaccard_idx_max else "↘️")
        print(f"Average similarity: {base_jaccard_idx_avg:.2f} {other_jaccard_idx_avg:.2f}", end=" ")
        print("↗️" if other_jaccard_idx_avg > base_jaccard_idx_avg else "↘️")

        other_oligosets = other_scores[region_id]
        compare_oligosets(base_oligosets, other_oligosets, relevant)


def main() -> None:
    base_oligosets_paths = read_oligosets_paths(BASE)
    other_oligosets_paths = read_oligosets_paths(OTHER)

    for pipeline_name, base_oligosets_path in base_oligosets_paths.items():
        if pipeline_name not in other_oligosets_paths:
            print(f"Skipping {pipeline_name}".upper())
            continue

        print(f"Comparing {pipeline_name}...".upper())
        relevant = "worst" if pipeline_name in WORST_RELEVANT else "average"
        compare_all_oligosets(base_oligosets_path, other_oligosets_paths[pipeline_name], relevant)
        print()


if __name__ == "__main__":
    main()
