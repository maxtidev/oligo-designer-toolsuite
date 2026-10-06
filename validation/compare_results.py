import os
from collections import defaultdict
from typing import Literal

import pandas as pd
import pyarrow as pa
import yaml

OUTPUT_PATH = "validation/outputs"

BASE = "baseline"
OTHER = "pyarrow"


def read_dfs_paths(
    subpath: str, variant: Literal["regions", "oligosets"] = "regions"
) -> dict[str, dict[str, list[str]]]:
    result = defaultdict(lambda: defaultdict(list))

    base_path = os.path.join(OUTPUT_PATH, subpath)
    for pipeline_output_dirname in os.listdir(base_path):
        pipeline_name = pipeline_output_dirname.split("_")[1]

        pipeline_output_path = os.path.join(base_path, pipeline_output_dirname)
        for database_dirname in os.listdir(pipeline_output_path):
            database_path = os.path.join(pipeline_output_path, database_dirname)
            if not os.path.isdir(database_path):
                continue

            steps_dirnames = sorted(
                [
                    dirname
                    for dirname in os.listdir(database_path)
                    if os.path.isdir(os.path.join(database_path, dirname))
                ]
            )

            result[pipeline_name][database_dirname] = [
                os.path.join(database_path, dirname, f"{variant}-snapshot.parquet")
                for dirname in steps_dirnames
            ]

    return dict(result)


def read_df(path: str) -> pd.DataFrame:
    df = pd.read_parquet(path, dtype_backend="pyarrow")
    # drop columns containing only null because of https://github.com/apache/arrow/issues/43838
    null_columns = [col for col in df if "null" in str(df[col].dtype) or df[col].count() == 0]
    df = df.drop(null_columns, axis=1)
    return df


def sorting_key_function(s: pd.Series) -> pd.Series:
    if "list" in str(s.dtype):
        # pandas/pyarrow doesn't support sorting by nested list, so we map to smallest value
        return pd.Series([min(min(value_list)) for value_list in s])
    return s


def wrap_series_in_list_of_lists(s: pd.Series) -> pd.Series:
    type = pa.array(s).type
    if type == pa.large_string():
        type = pa.string()

    nested_list_type = pa.list_(pa.field("element", pa.list_(pa.field("element", type))))
    new_values = [[[val if val is not pd.NA else None]] for val in s]
    return pd.Series(data=new_values, index=s.index, dtype=pd.ArrowDtype(nested_list_type))


def normalize_df(
    df: pd.DataFrame, sorting_keys: list[str], wrap_in_lists_of_lists: bool = False
) -> pd.DataFrame:
    sort_keys = []

    for key in sorting_keys:
        if key in df:
            sort_keys.append(key)

    assert sort_keys

    try:
        df = df.sort_values(sort_keys, kind="stable", key=sorting_key_function)
    except:
        print(df.dtypes)
        raise

    # drop index because oligo ids are arbitrary
    df = df.reset_index(drop=True)
    if "oligo_id" in df:
        df = df.drop("oligo_id", axis=1)

    # sort columns by name
    df = df.reindex(sorted(df.columns), axis=1)

    # ensure every value is in a list of lists
    if wrap_in_lists_of_lists:
        for column in df:
            if not "list" in str(df[column].dtype):
                df[column] = wrap_series_in_list_of_lists(df[column])
    return df


def assert_df_equality(base_df: pd.DataFrame, other_df: pd.DataFrame, step: str | None = None) -> bool:

    if other_df.equals(base_df):
        print(f"{step or ''} ", end="")
        return True
    else:
        print(f"{step or ''} MISMATCH ❌")
        for (base_name, base_col), (other_name, other_col) in zip(base_df.items(), other_df.items()):
            print(base_name, other_name)
            for (base_idx, base_val), (other_idx, other_val) in zip(base_col.items(), other_col.items()):
                if base_val is pd.NA and other_val is pd.NA:
                    continue
                if base_val != other_val:
                    print(
                        base_idx,
                        other_idx,
                        # pyrefly: ignore [bad-argument-type]
                        base_df["region_id"].iloc[int(base_idx)],
                        # pyrefly: ignore [bad-argument-type]
                        other_df["region_id"].iloc[int(other_idx)],
                        base_val,
                        other_val,
                    )
        print(len(base_df), len(other_df))

        return False


def assert_parquet_equality(base_parquet_path: str, other_parquet_path: str, step: str | None = None) -> bool:
    base_df = read_df(base_parquet_path)
    other_df = read_df(other_parquet_path)

    sorting_keys = ["region_id", "start", "oligo"]
    base_df = normalize_df(base_df, sorting_keys, wrap_in_lists_of_lists=True)
    other_df = normalize_df(other_df, sorting_keys, wrap_in_lists_of_lists=True)
    return assert_df_equality(base_df, other_df, step)


def read_output_paths(subpath: str) -> dict[str, list[str]]:
    result = defaultdict(list)

    base_path = os.path.join(OUTPUT_PATH, subpath)
    for pipeline_output_dirname in os.listdir(base_path):
        pipeline_name = pipeline_output_dirname.split("_")[1]

        pipeline_output_path = os.path.join(base_path, pipeline_output_dirname)
        for output_filename in os.listdir(pipeline_output_path):
            output_file_path = os.path.join(pipeline_output_path, output_filename)
            if os.path.isdir(output_filename):
                continue
            if output_filename.endswith("probes.yml"):
                # normalizing this one is too much of a hassle, rely on tsv comparison
                continue
            if any(output_filename.endswith(suffix) for suffix in (".yml", ".yaml", ".tsv", ".xlsx")):
                result[pipeline_name].append(output_file_path)

    return dict(result)


def assert_dict_equality(base_dict: dict, other_dict: dict) -> bool:
    if set(base_dict.keys()) != set(other_dict.keys()):
        return False
    for key, value in base_dict.items():
        other_value = other_dict[key]
        if isinstance(value, dict):
            if not isinstance(other_value, dict):
                print("strucural inequality")
                return False
            if not assert_dict_equality(value, other_value):
                return False
        else:
            if isinstance(other_value, dict):
                print("strucural inequality")
                return False
            if value != other_value:
                print(f"{key=}, {value=}, {other_value=}")
                return False
    return True


def assert_yml_equality(base_yml_path: str, other_yml_path: str) -> bool:
    # read into dict, sort keys and value dicts, compare
    with open(base_yml_path) as f:
        base_dict: dict = yaml.safe_load(f)
    with open(other_yml_path) as f:
        other_dict: dict = yaml.safe_load(f)

    return assert_dict_equality(base_dict, other_dict)


def assert_tsv_equality(base_tsv_path: str, other_tsv_path: str) -> bool:
    # read into df, sort, compare
    base_df = pd.read_csv(base_tsv_path, sep="\t")
    other_df = pd.read_csv(other_tsv_path, sep="\t")

    sorting_keys = ["region_id", "oligoset_id", "sequence_target", "target", "sequence_padlock_probe", "bit"]
    base_df = normalize_df(base_df, sorting_keys)
    other_df = normalize_df(other_df, sorting_keys)

    return assert_df_equality(base_df, other_df)


def assert_xlsx_equality(base_xlsx_path: str, other_xlsx_path: str) -> bool:
    # for each page: read into df, sort, compare
    return True


def assert_results_equality(base_output_path: str, other_output_path: str) -> bool:
    # check all output formats
    suffix = base_output_path.rsplit(".", maxsplit=1)[1]
    if suffix in ("yml", "yaml"):
        return assert_yml_equality(base_output_path, other_output_path)
    elif suffix == "tsv":
        return assert_tsv_equality(base_output_path, other_output_path)
    elif suffix == "xlsx":
        return assert_xlsx_equality(base_output_path, other_output_path)

    return False


def main():
    base_dfs_paths = read_dfs_paths(BASE)
    other_dfs_paths = read_dfs_paths(OTHER)

    base_output_paths = read_output_paths(BASE)
    other_output_paths = read_output_paths(OTHER)

    for pipeline_name, base_database_paths in base_dfs_paths.items():
        if pipeline_name not in other_dfs_paths:
            print(f"Skipping {pipeline_name}".upper())
            continue

        print(f"Comparing {pipeline_name}...".upper())
        for base_database_name, base_parquet_paths in base_database_paths.items():
            print(f"Comparing {base_database_name}...")
            other_parquet_paths = other_dfs_paths[pipeline_name][base_database_name]
            is_all_equal = True
            for i, (base_parquet_path, other_parquet_path) in enumerate(
                zip(base_parquet_paths, other_parquet_paths)
            ):
                is_equal = assert_parquet_equality(base_parquet_path, other_parquet_path, step=str(i + 1))
                is_all_equal = is_all_equal and is_equal
            if is_all_equal:
                print("EQUAL ✔️")
            else:
                print("NOT EQUAL ❌")

        for base_output_path, other_output_path in zip(
            base_output_paths[pipeline_name], other_output_paths[pipeline_name]
        ):
            print(f"Comparing {os.path.basename(base_output_path)}...")
            if assert_results_equality(base_output_path, other_output_path):
                print("EQUAL ✔️")
            else:
                print("NOT EQUAL ❌")


if __name__ == "__main__":
    main()
