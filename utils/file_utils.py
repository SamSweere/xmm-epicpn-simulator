import gzip
import os
import shutil
import tarfile

from utils import log


def check_file_exist(file_path):
    if os.path.exists(file_path):
        log.plog(f"File {file_path} already exists")
        return True
    else:
        return False


# Level 5 is a good balance between speed and compression ratio for my usecase
def compress_gzip(in_file_path, out_file_path, compresslevel=5):
    with open(in_file_path, 'rb') as f_in:
        with gzip.open(out_file_path, 'wb', compresslevel=compresslevel) as f_out:
            shutil.copyfileobj(f_in, f_out)


def decompress_gzip(in_file_path, out_file_dir):
    out_file_path = None
    with gzip.open(in_file_path, 'rb') as f_in:
        # out_file_path = os.path.join(out_file_dir
        out_file_path = os.path.join(out_file_dir, os.path.basename(f_in.filename)[:-3])  # Remove the .gz part
        with open(out_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return out_file_path


def compress_targz(in_file_path, out_file_path):
    with tarfile.open(out_file_path, "w:gz") as tar:
        tar.add(in_file_path, arcname=os.path.basename(in_file_path))


def decompress_targz(in_file_path, out_file_dir):
    decompressed_file_paths = []
    with tarfile.open(in_file_path, "r:gz") as tar:
        def is_within_directory(directory, target):
            
            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)
        
            prefix = os.path.commonprefix([abs_directory, abs_target])
            
            return prefix == abs_directory
        
        def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
        
            for member in tar.getmembers():
                member_path = os.path.join(path, member.name)
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")
        
            tar.extractall(path, members, numeric_owner=numeric_owner) 
            
        
        safe_extract(tar, path=out_file_dir)
        for i in range(len(tar.members)):
            name = tar.members[i].name
            decompressed_file_paths.append(os.path.join(out_file_dir, name))

    return decompressed_file_paths
