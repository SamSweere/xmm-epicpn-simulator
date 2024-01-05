import gzip
import os
import shutil
import tarfile
from pathlib import Path


# Level 5 is a good balance between speed and compression ratio for my usecase
def compress_gzip(in_file_path, out_file_path, compresslevel=5):
    with open(in_file_path, "rb") as f_in:
        with gzip.open(out_file_path, "wb", compresslevel=compresslevel) as f_out:
            shutil.copyfileobj(f_in, f_out)


def compress_targz(in_file_path: Path, out_file_path: Path):
    with tarfile.open(out_file_path, "w:gz") as tar:
        tar.add(in_file_path)


def decompress_targz(in_file_path: Path, out_file_dir: Path):
    out_file_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(in_file_path, "r:gz") as tar:

        def is_within_directory(directory: Path, target: Path):
            abs_directory = directory.absolute()
            abs_target = target.absolute()

            prefix = Path(os.path.commonpath([abs_directory, abs_target]))

            return prefix == abs_directory

        def safe_extract(tar: tarfile.TarFile, path: Path = Path(".")):
            for member in tar.getmembers():
                member_path = path / member.name
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")

            tar.extractall(path)

        safe_extract(tar, path=out_file_dir)
