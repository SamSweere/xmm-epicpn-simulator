from multiprocessing import cpu_count
from pathlib import Path
from typing import Annotated

from pydantic import (
    BaseModel,
    Field,
    NonNegativeFloat,
    NonNegativeInt,
    PositiveFloat,
    PositiveInt,
    computed_field,
)
from pydantic.functional_validators import AfterValidator


def _expanduser(v: Path) -> Path:
    return v.expanduser()


def _mkdir(v: Path) -> Path:
    v.mkdir(parents=True, exist_ok=True)
    return v


def _get_num_processes(v: int) -> int:
    if v == 0:
        return cpu_count()
    return v


def _check_xmm_filter(v: str) -> str:
    filters = ["thin", "med", "thick"]
    if v not in filters:
        raise ValueError(f"Unknown filter [{v}] for XMM-Newton! Available filters: {list(filters.keys())}")
    return v


def _check_modes(v: dict[str, list]) -> dict[str, list]:
    available_modes = ("proj", "slice")
    for mode in v:
        if mode not in available_modes:
            raise ValueError(f"Unknown mode '{mode}'! Available modes: {available_modes}")

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
ProcessCount = Annotated[NonNegativeInt, AfterValidator(_get_num_processes)]


class DownloadCfg(BaseModel):
    num_processes: ProcessCount
    top_n: PositiveInt
    resolutions: list[PositiveInt]
    snapshots: dict[Annotated[NonNegativeInt, Field(le=99)], PositiveFloat]
    simulations: dict[str, list[tuple[PositiveFloat, str]]]
    modes: Annotated[
        dict[str, list],
        AfterValidator(_check_modes),
    ]
    cutouts_path: CfgPath
    cutouts_compressed: Path
    fits_path: CfgPath
    fits_compressed: Path


# class _SimputModes(BaseModel):
#     img: Annotated[int, Field(ge=-1)]
#     agn: NonNegativeInt
#     bkg: bool


class _SimulationModes(BaseModel):
    img: Annotated[int, Field(ge=-1)]
    agn: Annotated[int, Field(ge=-1)]
    bkg: NonNegativeInt


class _SimputImg(BaseModel):
    n_gen: Annotated[int, Field(ge=-1)]


class _SimputAgn(BaseModel):
    n_gen: NonNegativeInt
    deblending_n_gen: NonNegativeFloat
    deblending_n_flux:NonNegativeFloat
    deblending_min_sep: NonNegativeFloat
    deblending_max_sep: NonNegativeFloat
    deblending_max_flux_delta: NonNegativeFloat
    put_source_in_center: NonNegativeInt


class _SimputBkg(BaseModel):
    n_gen: NonNegativeInt


class SimputCfg(BaseModel):
    num_processes: NonNegativeInt
    zoom_range: tuple[PositiveInt, PositiveInt]
    sigma_b_range: tuple[PositiveInt, PositiveInt]
    img: _SimputImg
    agn: _SimputAgn
    bkg: _SimputBkg
    # agn: dict[str, any]
    # bkg: dict[str, int]
    offset_std: PositiveFloat
    num_img_sample: PositiveInt
    simput_dir: CfgPath
    fits_dir: CfgPath
    fits_compressed: Path


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
    res_mults: list[PositiveInt]
    max_exposure: PositiveInt
    modes: _SimulationModes
    wait_time: NonNegativeFloat = 23.04e-6
    simput_dir: Path
    out_dir: CfgPath


class _EMOS(BaseModel):
    use: bool
    filter: Annotated[str, AfterValidator(_check_xmm_filter)]
    sim_separate_ccds: bool
    max_event_pattern: Annotated[int, Field(ge=-1, le=12)]

    @computed_field
    @property
    def filter_abbrv(self) -> str:
        filters = {"thin": "t", "med": "m", "thick": "k"}
        return filters[self.filter]


class _EPN(BaseModel):
    use: bool
    filter: Annotated[str, AfterValidator(_check_xmm_filter)]
    sim_separate_ccds: bool
    max_event_pattern: Annotated[int, Field(ge=-1, le=4)]

    @computed_field
    @property
    def filter_abbrv(self) -> str:
        filters = {"thin": "t", "med": "m", "thick": "k"}
        return filters[self.filter]


class XMMInstrument(BaseModel):
    emos1: _EMOS
    emos2: _EMOS
    epn: _EPN
