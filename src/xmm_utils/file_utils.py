import gzip
import shutil
from pathlib import Path
from src.xmm_utils.external_run import run_command


# Level 5 is a good balance between speed and compression ratio for my usecase
def compress_gzip(in_file_path, out_file_path, compresslevel=6):
    with open(in_file_path, "rb") as f_in:
        with gzip.open(out_file_path, "wb", compresslevel=compresslevel) as f_out:
            shutil.copyfileobj(f_in, f_out)


def compress_targz(in_file_path: Path, out_file_path: Path):
    run_command(f"tar -czf {out_file_path.resolve()} {in_file_path.resolve()}")


def decompress_targz(in_file_path: Path, out_file_dir: Path):
    out_file_dir.mkdir(parents=True, exist_ok=True)
    run_command(
        f"tar -xzf {in_file_path.resolve()} -C {out_file_dir.resolve()} --strip-components=1"
    )
