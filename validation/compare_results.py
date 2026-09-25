import os
import sys
from collections import defaultdict
from typing import Literal

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

OUTPUT_PATH = "validation/outputs"

BASE = "baseline"
OTHER = "pyarrow"


def read_paths(
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
    null_columns = [col for col in df if "null" in str(df[col].dtype)]
    df = df.drop(null_columns, axis=1)
    return df


def sorting_key_function(s: pd.Series) -> pd.Series:
    if "list" in str(s.dtype):
        # pandas/pyarrow doesn't support sorting by nested list, so we map to smallest value
        return pd.Series([min(min(value_list)) for value_list in s])
    return s


def wrap_series_in_list_of_lists(s: pd.Series) -> pd.Series:
    new_values = pa.array([[[val]] for val in s])
    return pd.Series(data=new_values, index=s.index, dtype=pd.ArrowDtype(new_values.type))


def normalize_df(df: pd.DataFrame) -> pd.DataFrame:
    if "start" in df:
        sort_keys = ["region_id", "start"]
    else:
        assert "oligo" in df
        sort_keys = ["oligo"]
    
    try:
        df = df.sort_values(sort_keys, kind="stable", key=sorting_key_function)
    except:
        print(df.dtypes)
        raise

    # drop index because oligo ids are arbitrary
    df = df.reset_index(drop=True)

    # sort columns by name
    df = df.reindex(sorted(df.columns), axis=1)

    # ensure every value is in a list of lists
    for column in df:
        if not "list" in str(df[column].dtype):
            df[column] = wrap_series_in_list_of_lists(df[column])
    return df


def assert_df_equality(base_parquet_path: str, other_parquet_path: str, step: str | None = None) -> bool:
    base_df = read_df(base_parquet_path)
    other_df = read_df(other_parquet_path)

    base_df = normalize_df(base_df)
    other_df = normalize_df(other_df)

    if other_df.equals(base_df):
        print(f"{step or ''} ", end="")
        return True
    else:
        # print(f"{step or ''} MISMATCH ❌")
        # return False
        print(base_df)
        print(other_df)

        print(base_df.dtypes)
        print(other_df.dtypes)
        sys.exit(1)


def main():
    base_dfs_paths = read_paths(BASE)
    other_dfs_paths = read_paths(OTHER)

    for pipeline_name, base_database_paths in base_dfs_paths.items():
        print(f"Comparing {pipeline_name}...".upper())
        for base_database_name, base_parquet_paths in base_database_paths.items():
            print(f"Comparing {base_database_name}...")
            other_parquet_paths = other_dfs_paths[pipeline_name][base_database_name]
            is_all_equal = True
            for i, (base_parquet_path, other_parquet_path) in enumerate(
                zip(base_parquet_paths, other_parquet_paths)
            ):
                is_equal = assert_df_equality(base_parquet_path, other_parquet_path, step=str(i + 1))
                is_all_equal = is_all_equal and is_equal
            if is_all_equal:
                print("EQUAL ✔️")
            else:
                print("NOT EQUAL ❌")


if __name__ == "__main__":
    main()
