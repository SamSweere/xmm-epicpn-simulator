import gzip
import os
import shutil
from pathlib import Path

from loguru import logger

from src.tools.cli import run_command


def compress_gzip(in_file_path: Path, out_file_path: Path, compresslevel=6, remove_file: bool = False):
    out_file_path.parent.mkdir(exist_ok=True, parents=True)
    with open(in_file_path, "rb") as f_in, gzip.open(out_file_path, "wb", compresslevel=compresslevel) as f_out:
        shutil.copyfileobj(f_in, f_out)

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
    if remove_files:
        shutil.rmtree(in_path)


def decompress_targz(in_file_path: Path, out_file_dir: Path, tar_options: str = ""):
    out_file_dir.mkdir(parents=True, exist_ok=True)
    run_command(f"tar -xzf {in_file_path.resolve()} -C {out_file_dir.resolve()} {tar_options}")
    logger.success(f"Decompressed {in_file_path} to {out_file_dir}")
