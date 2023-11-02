import os
from pathlib import Path
from typing import Literal


def get_ccf_path() -> Path:
    ccf_path = os.environ["SAS_CCFPATH"]
    if not ccf_path:
        raise NotADirectoryError(f"The environment variable 'SAS_CCFPATH' is not set!")

    ccf_path = Path(ccf_path)
    if not ccf_path.exists():
        raise NotADirectoryError(f"The CCF directory does not exist at {ccf_path.resolve()}")
    return ccf_path


def get_epn_lincoord() -> Path:
    ccf_path = get_ccf_path()
    epn_lincoord = list(ccf_path.glob("EPN_LINCOORD*.CCF"))
    if not epn_lincoord:
        raise FileExistsError(f"Could not find any EPN_LINCOORD*.CCF in '{ccf_path.resolve()}'!")
    epn_lincoord.sort()
    return epn_lincoord[-1]


def get_emos_lincoord(emos_num: Literal[1, 2]):
    ccf_path = get_ccf_path()
    emos_lincoord = list(ccf_path.glob(f"EMOS{emos_num}_LINCOORD*.CCF"))
    if not emos_lincoord:
        raise FileExistsError(f"Could not find any EMOS{emos_lincoord}_LINCOORD*.CCF in '{ccf_path.resolve()}'!")
    emos_lincoord.sort()
    return emos_lincoord[-1]


def get_xmm_miscdata() -> Path:
    ccf_path = get_ccf_path()
    xmm_miscdata = list(ccf_path.glob("XMM_MISCDATA*.CCF"))
    if not xmm_miscdata:
        raise FileExistsError(f"Could not find any XMM_MISCDATA*.CCF in '{ccf_path.resolve()}'!")
    xmm_miscdata.sort()
    return xmm_miscdata[-1]
