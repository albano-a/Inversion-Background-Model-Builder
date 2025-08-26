import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import segyio
import gstools as gs
import os
import lasio
from pathlib import Path
import re
from scipy.interpolate import interp1d, interp2d
from scipy.interpolate import griddata, RegularGridInterpolator, LinearNDInterpolator
from scipy.ndimage import map_coordinates
from scipy import signal
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline

FT_TO_US = 304800


def open_segyio_cube(filepath, iline=189, xline=193):
    seis = segyio.open(filepath, iline=iline, xline=xline)
    vel = segyio.cube(seis)
    vel = vel[:-1, :-1, 425:1025]  # Hardcoded
    return vel


def create_welldata_matrix(folder: str) -> pd.DataFrame:
    """This function generates the data matrix for the wells

    Args:
        folder (str): path to the folder with well logs

    Returns:
        pd.DataFrame: Dataframe with (nwells, nsamples)
    """
    pattern = re.compile(r"~\s*DEPTH", re.IGNORECASE)
    well_skiprows, well_data = {}, {}

    for file in Path(folder).glob("*.las"):
        well_name = file.stem
        well_skiprows[well_name] = 0
        with file.open("r", encoding="utf-8", errors="ignore") as f:
            for i, line in enumerate(f):
                if pattern.match(line):
                    well_skiprows[well_name] = i + 1
                    break

    for file in Path(folder).glob("*.las"):
        well_name = file.stem
        skiprows = well_skiprows[well_name]

        names = pd.read_csv(
            file,
            sep=r"\s+",
            skiprows=skiprows - 1,
            nrows=0,
        )
        names = names.columns[1:].str.strip()

        df = pd.read_csv(
            file,
            sep=r"\s+",
            skiprows=skiprows,
            na_values=-9999,
            usecols=lambda c: c in names,
            names=names,
        )

        df["CALIR"] = np.nan

        if "DTCO" in df:
            df["VP"] = FT_TO_US / df["DTCO"]
            df["IP"] = df["VP"] * df["RHOB"]
        else:
            df["VP"], df["IP"] = np.nan, np.nan

        if "DTS" in df:
            df["VS"] = FT_TO_US / df["DTS"]
            df["IS"] = df["VS"] * df["RHOB"]
        else:
            df["VS"], df["IS"] = np.nan, np.nan

        df["Well"] = well_name
        well_data[well_name] = df

    return pd.concat(well_data.values(), ignore_index=True)


def create_wellposition_matrix(dataframe, time):
    """
    Considere uma matriz com os poços empilhados um por linha contendo o mesmo tamanho e amostrado
    na mesma dimensão do modelo,
    ou seja: (nwells, nsamples)
    """
    log_well = np.zeros(shape=(16, 600))  # Hardcoded
    i = 0
    dataframe["IP_ups"] = dataframe["IP"]
    for well in dataframe.Well.unique():
        mask = dataframe.Well == well
        dataframe.loc[mask, "IP_ups"] = np.convolve(
            np.ones(30) / 30, dataframe[mask].IP, mode="same"
        )
        interpolator = interp1d(
            dataframe[mask].TWT, dataframe[mask].IP_ups, bounds_error=False, fill_value='extrapolate'
        )
        log_well[i,:] = interpolator(time)
        i += 1


def find_normalized_thickness(surfs, loc):
    """
    Compute the normalized cumulative thickness between surfaces.
    Parameters
    ----------
    surfs : list of ndarray
        List of 2D arrays representing surfaces.
    loc : tuple of int
        Location in inline and crossline indices (i_il, i_xl).
    Returns
    -------
    ndarray
        Cumulative normalized thickness.
    """
    i_il, i_xl = loc
    thick = np.zeros(len(surfs))
    for i in range(1, len(surfs)):
        thick[i] = abs(surfs[i][i_il, i_xl] - surfs[i - 1][i_il, i_xl])
    thick = abs(thick) / np.sum(thick)
    return np.cumsum(thick)


def deflat(cube, ijk, silence=True):
    """
    Deflat a data cube using precomputed indices.
    Parameters
    ----------
    cube : ndarray
        3D seismic data array (il,xl,nsamples).
    ijk : ndarray
        3D index map for deflattening (il,xl,nsamples).
    silence : bool, optional
        If True, suppress progress messages. Default is True.
    Returns
    -------
    ndarray
        Deflattened cube.
    """
    new_cube = np.empty_like(ijk)
    for i in range(cube.shape[0]):
        cube_cut = cube[i].T.copy()
        ijk_cut = ijk[i].T.copy()
        _, x = np.indices(ijk_cut.shape)
        y = ijk_cut * (cube.shape[2] - 1)
        y[:, -1] = y[:, -2]
        deformed = map_coordinates(cube_cut, (y, x))
        new_cube[i] = deformed.T
        if not silence:
            print(
                f"Deflattening Inline Index [{i}] {((i + 1)/ cube.shape[0]) * 100:.2f}%\r",
                end="",
            )
    return new_cube


def low_pass(cube, cut_hz_vert, dt_s, cut_hz_hor=15, silence=True):
    """
    Apply a low-pass filter to a data cube.
    Parameters
    ----------
    cube : ndarray
        3D seismic data array (il,xl,nsamples).
    cut_hz : float
        Cut-off frequency in Hz.
    dt_s : float
        Sampling interval in seconds.
    silence : bool, optional
        If True, suppress progress messages. Default is True.
    Returns
    -------
    ndarray
        Filtered seismic cube.
    """
    nfilt = 1
    samplerate_hz = 1 / dt_s
    nyquist = samplerate_hz / 2
    b, a = signal.butter(nfilt, cut_hz_vert / nyquist, btype="lowpass")
    b2, a2 = signal.butter(nfilt, cut_hz_hor / nyquist, btype="lowpass")
    new_cube = np.empty_like(cube)
    for i in range(cube.shape[0]):
        new_cube[i, :, :] = signal.filtfilt(b, a, cube[i, :, :], axis=0)
        new_cube[i, :, :] = signal.filtfilt(b2, a2, cube[i, :, :], axis=1)
        if not silence:
            print(
                f"Filtering Inline Index [{i}] {((i + 1)/ cube.shape[0]) * 100:.2f}%\r",
                end="",
            )
    for j in range(cube.shape[1]):
        new_cube[:, j, :] = signal.filtfilt(b2, a2, cube[:, j, :], axis=1)
        if not silence:
            print(
                f"Filtering Crossline Index [{j}] {((j + 1)/ cube.shape[1]) * 100:.2f}%\r",
                end="",
            )
    return new_cube


class BMTools:
    """
    Class for handling a full background modeling workflow.
    """

    def __init__(self):
        """Initialize BMTools with default attributes."""
        self.logs = None  # Well logs
        self.pos = None  # Relative well locations
        self.seismic_shape = None  # Shape of the seismic cube
        self.rgt_thick = None  # Relative Geological Time thickness
        self.rgt = None  # RGT cube
        self.logs_flat = None  # Flattened logs
        self.krig_flat = None  # Kriging results
        self.krig = None  # Deformed cube
        self.time = None  # Time vector
        self.dt = None  # Time sampling interval
        self.dt_s = None  # Time sampling in seconds
        self.rgt_ref = None  # Reference RGT trace
        self.model = None  # Low-pass filtered data (background model)

    def build_rgt(self, surfaces, time, silence=True, ref_trace=None):
        """
        Build Relative Geological Time (RGT) volume.
        Parameters
        ----------
        surfaces : ndarray
            3D array of interpreted surfaces (il,xl,nsamples).
        time : ndarray
            1D array of time values in milliseconds (nsamples).
        silence : bool, optional
            If True, suppress progress messages. Default is True.
        Returns
        -------
        BMTools
            Self instance with updated RGT attributes.
        """
        self.seismic_shape = np.asarray(
            (surfaces.shape[1], surfaces.shape[2], len(time))
        )
        if not ref_trace:
            ref_trace = (self.seismic_shape[0] // 2, self.seismic_shape[1] // 2)
        self.rgt_thick = find_normalized_thickness(surfaces, ref_trace)
        self.rgt = np.full(self.seismic_shape, np.nan)
        self.time = time.copy()
        self.dt = self.time[1] - self.time[0]
        self.dt_s = self.dt / 1000
        nlayers = len(surfaces) - 1

        for i in range(self.seismic_shape[0]):
            for j in range(self.seismic_shape[1]):
                for k in range(nlayers):
                    mask = np.arange(
                        (surfaces[k][i, j] - self.time[0]) // self.dt,
                        (surfaces[k + 1][i, j] - self.time[0]) // self.dt,
                    ).astype(int)
                    self.rgt[i, j, mask] = np.linspace(
                        self.rgt_thick[k], self.rgt_thick[k + 1], len(mask)
                    )
                self.rgt[i, j, : np.nanargmin(self.rgt[i, j])] = 0
                self.rgt[i, j, np.nanargmax(self.rgt[i, j]) :] = 1
            if not silence:
                print(
                    f"Building RGT Model Inline Index [{i}] ({((i + 1)/ self.seismic_shape[0]) * 100:.2f})%\r",
                    end="",
                )
        return self

    def read_wells(self, well_logs, relative_well_loc, rgt_domain_length=500):
        """
        Process well logs and map to RGT.
        Parameters
        ----------
        well_logs : ndarray
            2D array of well logs (nwells, nsamples).
        relative_well_loc : ndarray
            2D array of well locations (index of inline, index of crossline).
        Returns
        -------
        BMTools
            Self instance with updated well log attributes.
        """
        self.logs_z = well_logs
        self.pos = relative_well_loc.T
        self.logs_flat = np.empty((self.logs_z.shape[0], rgt_domain_length))
        self.rgt_ref = np.linspace(0, 1, rgt_domain_length)

        for i in range(self.logs_flat.shape[0]):
            il_index, xl_index = self.pos
            index_rgt = self.rgt[il_index[i], xl_index[i], :].copy()
            index_rgt[: np.nanargmin(index_rgt)] = 0
            index_rgt[np.nanargmax(index_rgt) + 1 :] = 1
            interpolator = interp1d(
                index_rgt, self.logs_z[i], fill_value=np.nan, kind="linear"
            )
            self.logs_flat[i, :] = interpolator(self.rgt_ref)
        return self

    def interpol(
        self,
        variogram,
        well_logs,
        relative_well_loc,
        fill_value,
        skip_values=0,
        second_var=None,
        silence=True,
        decimate=0,
        rgt_domain_length=500,
    ):
        """
        Perform kriging, deflattening and low-pass filter.
        Parameters
        ----------
        variogram : gstools.Variogram
            Variogram (gsstools) model for kriging.
        cut_hz : float
            Cut-off frequency in Hz for low-pass filtering.
        fill_value : tuple of float
            Values to fill in top and base regions.
        silence : bool, optional
            If True, suppress progress messages. Default is True.
        Returns
        -------
        BMTools
            Self instance with updated attributes.
        """
        self.read_wells(well_logs, relative_well_loc, rgt_domain_length)
        self.rgt_ref = np.linspace(0, 1, rgt_domain_length)
        self.krig_flat = np.zeros(
            (self.rgt.shape[0], self.rgt.shape[1], rgt_domain_length)
        )
        self.krig = np.zeros_like(self.rgt)
        gridy = np.arange(self.seismic_shape[0], dtype=np.float32)
        gridx = np.arange(self.seismic_shape[1], dtype=np.float32)
        x, y = (self.pos[1], self.pos[0])
        top, base = True, False
        if decimate > 0:
            gridx_dec = gridx.copy()[::decimate]
            gridy_dec = gridy.copy()[::decimate]

        for i in range(self.rgt_ref.shape[0]):
            if self.rgt_ref[i] == 0:
                self.krig_flat[:, :, i] = fill_value[0]
            elif self.rgt_ref[i] == 1:
                self.krig_flat[:, :, i] = fill_value[1]
            else:
                dvals = self.logs_flat[:, i].copy()
                # dvals[dvals==skip_values] = np.nan
                if np.sum(~np.isnan(dvals)) > 0:
                    OK = gs.krige.ExtDrift(
                        variogram,
                        (x, y),
                        dvals,
                        ext_drift=second_var[
                            self.pos[0][~np.isnan(dvals)],
                            self.pos[1][~np.isnan(dvals)],
                            i,
                        ],
                    )
                    if decimate <= 0:
                        self.krig_flat[:, :, i] = OK.structured(
                            [gridx, gridy], ext_drift=second_var[:, :, i]
                        )[0].T
                    else:
                        krig_dec = OK.structured(
                            [gridx_dec, gridy_dec],
                            ext_drift=second_var[::decimate, ::decimate, i],
                        )[0]
                        interpol = RectBivariateSpline(gridx_dec, gridy_dec, krig_dec)
                        OK_structured = interpol(gridx, gridy)
                        self.krig_flat[:, :, i] = OK_structured.T
                    top, base = False, True
                elif top:
                    self.krig_flat[:, :, i] = fill_value[0]
                elif base:
                    self.krig_flat[:, :, i] = fill_value[1]
            if not silence:
                print(
                    f"Kriging Index [{i}] ({((i + 1) / self.seismic_shape[2]) * 100:.2f}%) done\r",
                    end="",
                )
        self.krig = deflat(self.krig_flat, self.rgt, silence)
        return self

    def smooth_model(self, cut_hz_vert, cut_hz_hor=15, silence=True):
        self.model = low_pass(self.krig, cut_hz_vert, self.dt_s, cut_hz_hor, silence)
        return self
