import os
from pathlib import Path

from src.xmm_utils.external_run import run_command


def compress_gzip(in_file_path: Path, out_file_path: Path, compresslevel=6, remove_file: bool = False):
    run_command(f"gzip -{compresslevel} -c {in_file_path.resolve()} > {out_file_path.resolve()}")
    if remove_file:
        in_file_path.unlink()


def compress_targz(in_path: Path, out_file_path: Path, remove_files: bool = False):
    if not out_file_path.name.endswith(".tar.gz"):
        raise ValueError(f"Output file path {out_file_path.resolve()} does not end with '.tar.gz'")
    out_file_path.parent.mkdir(parents=True, exist_ok=True)
    suffix = " --remove-files" if remove_files else ""
    run_command(
        f"cd {in_path.parent.resolve()} && "
        + f"tar -czf {out_file_path.resolve()} {in_path.name}{os.sep} --overwrite{suffix}"
    )


def decompress_targz(in_file_path: Path, out_file_dir: Path):
    out_file_dir.mkdir(parents=True, exist_ok=True)
    run_command(f"tar -xzf {in_file_path.resolve()} -C {out_file_dir.resolve()} --strip-components=1")
