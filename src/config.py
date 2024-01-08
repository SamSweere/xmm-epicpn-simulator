from typing import Annotated, Dict, List, Literal, Tuple
from pathlib import Path
from pydantic import (
    BaseModel,
    Field,
    PositiveInt,
    NonNegativeFloat,
    NonNegativeInt,
    PositiveFloat,
)
from pydantic.functional_validators import AfterValidator
from multiprocessing import cpu_count


def _expanduser(v: Path) -> Path:
    return v.expanduser()


def _mkdir(v: Path) -> Path:
    v.mkdir(parents=True, exist_ok=True)
    return v


def _get_num_processes(v: int) -> int:
    if v == 0:
        return cpu_count()
    return v


def _check_modes(v: Dict[str, list]) -> Dict[str, list]:
    available_modes = ("proj", "slice")
    for mode in v.keys():
        if mode not in available_modes:
            raise ValueError(
                f"Unknown mode '{mode}'! Available modes: {available_modes}"
            )

        axes = v[mode]
        for axis in axes:
            if axis not in ("x", "y", "z"):
                raise ValueError(f"Unknown axis '{axis}'! Available axes: x, y, z")

        if mode == "proj":
            new_axes = []
            for axis in axes:
                if axis == "x":
                    new_axes.append([1, 0, 0])
                if axis == "y":
                    new_axes.append([0, 1, 0])
                if axis == "z":
                    new_axes.append([0, 0, 1])
            v[mode] = new_axes

    return v


CfgPath = Annotated[Path, AfterValidator(_expanduser), AfterValidator(_mkdir)]
InstrumentName = Literal["epn", "emos1", "emos2"]
XMMFilter = Literal["thin", "med", "thick"]
ProcessCount = Annotated[NonNegativeInt, AfterValidator(_get_num_processes)]


class DownloadCfg(BaseModel):
    num_processes: ProcessCount
    top_n: PositiveInt
    resolutions: List[PositiveInt]
    snapshots: Dict[Annotated[NonNegativeInt, Field(le=99)], PositiveFloat]
    simulations: Dict[str, List[Tuple[PositiveFloat, str]]]
    modes: Annotated[
        Dict[str, list],
        AfterValidator(_check_modes),
    ]
    cutouts_path: CfgPath
    cutouts_compressed: Path
    fits_path: CfgPath
    fits_compressed: Path


class _SimputModes(BaseModel):
    img: Annotated[int, Field(ge=-1)]
    agn: NonNegativeInt
    bkg: bool


class _SimulationModes(BaseModel):
    img: Annotated[int, Field(ge=-1)]
    agn: Annotated[int, Field(ge=-1)]
    bkg: NonNegativeInt


class SimputCfg(BaseModel):
    num_processes: NonNegativeInt
    instruments: List[InstrumentName]
    filter: XMMFilter
    zoom_range: Tuple[PositiveInt, PositiveInt]
    sigma_b_range: Tuple[PositiveInt, PositiveInt]
    modes: _SimputModes
    offset_std: PositiveFloat
    num_img_sample: PositiveInt
    simput_dir: CfgPath
    fits_dir: CfgPath
    fits_compressed: Path
    agn_dir: CfgPath
    bkg_dir: CfgPath
    img_dir: CfgPath


class EnergySettings(BaseModel):
    emin: NonNegativeFloat
    emax: NonNegativeFloat


class EnvironmentCfg(BaseModel):
    working_dir: CfgPath
    output_dir: CfgPath
    log_dir: CfgPath
    fail_on_error: bool = False
    debug: bool = False
    verbose: bool = True
    overwrite: bool = False
    consume_data: bool = False


class SimulationCfg(BaseModel):
    num_processes: NonNegativeInt
    instruments: List[InstrumentName]
    filter: XMMFilter
    res_mults: List[PositiveInt]
    max_exposure: PositiveInt
    modes: _SimulationModes
    sim_separate_ccds: bool
    wait_time: NonNegativeFloat = 23.04e-6
    simput_dir: Path
    out_dir: CfgPath