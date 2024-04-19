from pathlib import Path


def run(path_to_cfg: Path):
    with open(path_to_cfg) as file:  # noqa
        pass


if __name__ == "__main__":
    run(Path())
