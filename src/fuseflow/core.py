"""Core functionality for MarrmotFlow package."""
# internal imports
from ._default_dicts import (
    default_forcing_units,
    default_forcing_vars,
    default_models,
)

from .templating import render_settings

# built-in imports
import glob
import os
import warnings
import json
import shutil
import sys
import re

from typing import (
    Dict,
    Sequence,
    Union
)
from dateutil import parser


# third-party imports
import xarray as xr
import cdo
import geopandas as gpd
import pandas as pd
import pint_xarray
import pint
import pyet
import numpy as np
import timezonefinder

# define types
try:
    from os import PathLike
except ImportError:
    PathLike = str
else:
    PathLike = Union[str, PathLike]

# global variables
pkgdir = sys.modules['fuseflow'].__path__[0]
setting_path = os.path.join(pkgdir, 'default_settings')


class FUSEWorkflow:
    """A class representing a FUSE workflow."""

    def __init__(
        self,
        name: str = "fuse_models",
        cat: gpd.GeoDataFrame | PathLike = None, # type: ignore
        forcing_vars: Dict[str, str] = {},
        forcing_files: Sequence[PathLike] | PathLike = None, # type: ignore
        forcing_units: Dict[str, str] = {},
        pet_method: str = "hamon",
        model_number: Sequence[int] | int = default_models, # HBV-96 and GR4J as default models
        forcing_time_zone: str = None,
        model_time_zone: str = None,
        streamflow: xr.DataArray | PathLike = None, # type: ignore
        elev_bands: PathLike | str = None, # type: ignore
        settings: Dict = {},
    ) -> 'FUSEWorkflow':
        """
        Initialize the FUSE workflow with forcing variables, files, units, PET method, and model number.

        Parameters
        ----------
        forcing_vars : Dict[str, str], optional
            Dictionary of forcing variables and their units, by default {}.
            The mandatory keys to this dictionary are:
            - 'precip': Precipitation variable name
            - 'temp': Temperature variable name
            The values must be present in the `forcing_files`.
        forcing_files : Sequence[PathLike] | PathLike, optional
            Sequence of paths to forcing files or a single path, by default None
        forcing_units : Dict[str, str], optional
            Dictionary of units for the forcing variables, by default {}.
            The keys must match the keys in `forcing_vars`, and the values
            must be valid pint units.
        pet_method : str, optional
            Method for calculating potential evapotranspiration, by default "hamon"
        model_number : Sequence[int] | int, optional
            Model number(s) to be used in the workflow, by default [7, 37] (HBV-96 and GR4J)
        forcing_time_zone : str, optional
            Time zone for the forcing data, by default None.
            If None, it will default to 'UTC'.
        model_time_zone : str, optional
            Time zone for the model, by default None.
            If None, it will be autodetected based on the catchment's centroid coordinates.
        cat : gpd.GeoDataFrame | PathLike, optional
            Catchment data as a GeoDataFrame or a path to a GeoJSON file, by default None.
            If None, it raises a ValueError.
        streamflow : xr.DataArray | PathLike, optional
            Streamflow data as a xarray.DataArray, by default None. For now, only the
            NetCDF format in form of a xarray.DataArray is supported. The observation
            must be in the 'q_obs' variable of the NetCDF file.

        Raises
        ------
        ValueError
            If forcing files are not provided or are not in the correct format.
        TypeError
            If forcing files are not a sequence or a PathLike object.

        Notes
        -----
        - `pet_method` only accepts "hamon" as a valid method. Other methods
           will be added in the future.
        """
        # assing the name of the workflow
        self.name = name

        # assign the catchment (cat) as a GeoDataFrame or PathLike
        if cat is None:
            raise ValueError("Catchment (cat) must be provided as a GeoDataFrame or PathLike.")
        if isinstance(cat, gpd.GeoDataFrame):
            self.cat = cat
        elif isinstance(cat, PathLike):
            self.cat = gpd.read_file(cat)
        else:
            raise TypeError("cat must be a GeoDataFrame or a PathLike object.")
        
        # assign forcing variables
        self.forcing_vars = forcing_vars

        # if not a list of forcing files and is a PathLike, read the files using glob.glob
        if isinstance(forcing_files, PathLike):
            forcing_files = glob.glob(os.path.join(forcing_files, "**/*.nc*"), recursive=True)

        # if not, then the user must provide a list of forcing files in NetCDF format
        # assigning a class attribute
        self.forcing_files = forcing_files

        # assign forcing units
        self.forcing_units = forcing_units

        # assign pet method
        self.pet_method = pet_method

        # assign model number
        if isinstance(model_number, int):
            model_number = [model_number]
        self.model_number = model_number

        # assign time zones
        self.forcing_time_zone = forcing_time_zone
        self.model_time_zone = model_time_zone

        # streamflow data
        if isinstance(streamflow, PathLike):
            # if streamflow is a PathLike, read it as a DataArray
            ds = xr.open_dataset(streamflow, engine='netcdf4')
            self.streamflow = ds['q_obs']
        elif isinstance(streamflow, xr.DataArray):
            # if streamflow is a DataArray, assign it directly
            self.streamflow = streamflow
        else:
            self.streamflow = None

        # settings for the workflow
        if not isinstance(settings, dict):
            raise TypeError("settings must be a dictionary.")
        self.settings = settings
        # check if mandatory values are provided
        mandatory_settings = ['start_date', 'end_date']
        for key in mandatory_settings:
            if key not in self.settings:
                raise ValueError(f"Missing mandatory setting: {key}")
            
        # assign elevation bands
        self.elev_bands = elev_bands

        # assign an output object
        self.output_mat = None  # Placeholder for output matrix

    @classmethod
    def from_dict(
        cls: 'FUSEWorkflow',
        init_dict: Dict = {},
    ) -> 'FUSEWorkflow':
        """
        Constructor to use a dictionary to instantiate
        """
        if len(init_dict) == 0:
            raise KeyError("`init_dict` cannot be empty")
        assert isinstance(init_dict, dict), "`init_dict` must be a `dict`"

        return cls(**init_dict)

    @classmethod
    def from_json(
        cls: 'FUSEWorkflow',
        json_str: str,
    ) -> 'FUSEWorkflow':
        """
        Constructor to use a loaded JSON string
        """
        # building customized FUSEWorkflow's JSON string decoder object
        decoder = json.JSONDecoder(object_hook=FUSEWorkflow._json_decoder)
        json_dict = decoder.decode(json_str)
        # return class instance
        return cls.from_dict(json_dict)

    @classmethod
    def from_json_file(
        cls: 'FUSEWorkflow',
        json_file: 'str',
    ) -> 'FUSEWorkflow':
        """
        Constructor to use a JSON file path
        """
        with open(json_file) as f:
            json_dict = json.load(f, object_hook=FUSEWorkflow._json_decoder)

        return cls.from_dict(json_dict)

    @staticmethod
    def _env_var_decoder(s):
        """
        OS environmental variable decoder
        """
        # RE patterns
        env_pat = r'\$(.*?)/'
        bef_pat = r'(.*?)\$.*?/?'
        aft_pat = r'\$.*?(/.*)'
        # strings after re matches
        e = re.search(env_pat, s).group(1)
        b = re.search(bef_pat, s).group(1)
        a = re.search(aft_pat, s).group(1)
        # extract environmental variable
        v = os.getenv(e)
        # return full: before+env_var+after
        if v:
            return b+v+a
        return s

    @staticmethod
    def _is_valid_integer(obj):
        """
        Check if a string can be converted to an integer
        """
        try:
            int(obj)
            return True
        except (ValueError, TypeError):
            return False

    @staticmethod
    def _json_decoder(obj):
        """
        Decoding typical JSON strings returned into valid Python objects
        """
        if obj in ["true", "True", "TRUE"]:
            return True
        elif obj in ["false", "False", "FALSE"]:
            return False
        elif isinstance(obj, str):
            if '$' in obj:
                return FUSEWorkflow._env_var_decoder(obj)
            if FUSEWorkflow._is_valid_integer(obj):
                return int(obj)
        elif isinstance(obj, dict):
            return {FUSEWorkflow._json_decoder(k): FUSEWorkflow._json_decoder(v) for k, v in obj.items()}
        return obj

    # class methods
    def run(self):
        """Run the workflow."""
        self.init_forcing_files() # defines self.df
        self.init_pet() # defines self.pet
        self.init_streamflow() # defines self.forcing['q_obs']
        self.init_elev_bands() # defines self.elev_bands

        # print a message about the timezones
        print(f"Using forcing time zone: {self.forcing_time_zone}")
        print(f"Using model time zone: {self.model_time_zone}")

        return f"Workflow executed successfully with {len(self.forcing_files)} forcing files."

    def save(self, save_path: PathLike): # type: ignore
        """Save the workflow output to a specified path."""
        if not hasattr(self, 'forcing') or self.forcing is None:
            raise ValueError("No output matrix to save. Run the workflow first.")
        
        # check if the save_path exists, if not, create it
        if not os.path.exists(save_path):
            os.makedirs(save_path, exist_ok=True)
        # also check for the save_path/settings, input, and output directories
        os.makedirs(os.path.join(save_path, 'settings'), exist_ok=True)
        os.makedirs(os.path.join(save_path, 'input'), exist_ok=True)
        os.makedirs(os.path.join(save_path, 'output'), exist_ok=True)

        # Create .mat file using the scipy.io.savemat function
        # the dataframe must be a cobination of self.df and self.pet
        for model in self.model_number:
            content = self.init_model_file(base_path=save_path, model_n=model)
            # save the `fm_catch` content to a text file
            with open(os.path.join(save_path, f'{self.name}_{model}.txt'), 'w') as f:
                f.write(content)

        # save the forcing data to a NetCDF file
        self.forcing.to_netcdf(os.path.join(save_path, 'input', f'{self.name}_input.nc'))
        # save the elevation bands to a NetCDF file
        self.elev_bands.to_netcdf(os.path.join(save_path, 'input', f'{self.name}_elev_bands.nc'))

        # copy the defaults files and folders to the settings directory
        for f in glob.glob(os.path.join(setting_path, '*')):
            # if a file, copy it to the settings directory
            if os.path.isfile(f):
                shutil.copy2(f, os.path.join(save_path, 'settings'))
            # if a directory, copy it to the settings directory
            elif os.path.isdir(f):
                shutil.copytree(f, os.path.join(save_path, 'settings', os.path.basename(f)))

        return f"Outputs saved to {save_path}"

    def __str__(self):
        return f"FUSEWorkflow(forcing_files={self.forcing_files})"
    
    def init_forcing_files(self):
        """Initialize forcing files."""
        # check the timezones
        if self.forcing_time_zone is None:
            warnings.warn(
                "Forcing time zone is not set. Defaulting to 'UTC'.",
                UserWarning
            )
            self.forcing_time_zone = "UTC"

        # check the model time-zone
        if self.model_time_zone is None:
            # try extracting the time zone from the provided `cat` file
            # calculate area's centroid coordinates---basically what is done
            # in the `init_class` method
            warnings.warn(
                "No `model_time_zone` provided in the settings. "
                "Autodetecting the time zone using `timezonefinder` "
                "based on the centroid coordinates of the catchment.",
                UserWarning,
            )
            if not self.cat.crs:
                warnings.warn(
                    "Catchment (cat) does not have a CRS. A default of "
                    "EPSG:4326 will be used for latitude calculation.",
                    UserWarning)
                self.cat.set_crs(epsg=4326, inplace=True)

            # if crs is not set to EPSG:4326, then convert it to EPSG:4326
            elif self.cat.crs != 'EPSG:4326':
                self.cat.to_crs(epsg=4326, inplace=True)

            # calculate the latitude from the catchment geometry
            lat = self.cat.geometry.centroid.y.mean()
            lng = self.cat.geometry.centroid.x.mean()

            # extracing the model time zone from the coordinates
            self.model_time_zone = timezonefinder.TimezoneFinder().timezone_at(
                lat=lat,
                lng=lng
            )

            # Print the model time zone
            if self.model_time_zone:
                warnings.warn(
                    f"Autodetected model time zone: {self.model_time_zone}",
                    UserWarning,
                )
            # if the model time zone is None, then assume UTC
            # and warn the user
            else:
                self.model_time_zone = 'UTC'
                warnings.warn(
                    "No `model_time_zone` provided in the settings and"
                    " autodetection using `timezonefinder` failed."
                    " Assuming UTC time zone.",
                    UserWarning,
                )

        _ureg = pint.UnitRegistry(force_ndarray_like=True)

        # read the forcing files using xarray and create a self.forcing
        cdo_obj = cdo.Cdo()  # CDO object
        ds = cdo_obj.mergetime(input=self.forcing_files, returnXArray=list(self.forcing_vars.values()))  # Mergeing

        # adjust the model time zone
        if self.model_time_zone != self.forcing_time_zone:
            # convert the self.forcing to the forcing time zone
            ds = ds.assign_coords({
                   'time': ds.time.to_index().tz_localize(self.forcing_time_zone).tz_convert(self.model_time_zone).tz_localize(None)
                })

        # rename the self.forcing variables to match the forcing_vars
        # and assign pint units to the self.forcing variables
        rename_vars_dict = {}
        for key, value in self.forcing_vars.items():
            if value in ds.variables:
                rename_vars_dict[value] = default_forcing_vars.get(key)
        ds = ds.rename(rename_vars_dict)

        # drop the variable not in self.forcing_vars
        ds = ds[list(default_forcing_vars.values())]

        # assign pint units to the self.forcing variables
        renamed_forcing_units = {}
        for key, value in self.forcing_units.items():
            new_key = default_forcing_vars.get(key)
            renamed_forcing_units[new_key] = value
        ds = ds.pint.quantify(units=renamed_forcing_units, unit_registry=_ureg)

        # convert the self.forcing units to the default forcing units
        renamed_to_forcing_units = {}
        for key, value in default_forcing_units.items():
            new_key = default_forcing_vars.get(key)
            renamed_to_forcing_units[new_key] = value
        ds = ds.pint.to(units=renamed_to_forcing_units)

        # after unit conversion, dequantify the self.forcing
        ds = ds.pint.dequantify()

        # calculate the centroid coordinates of the catchment(s)
        # Extract COMID and centroid coordinates
        # the id should be common between `self.cat` and `ds`
        id_to_coords = {}
        # `id` is basically the dimension not named `time` in `ds`
        id = list(set(ds.dims) - {'time'})[0]
        for idx, row in self.cat.iterrows():
            centroid = row['geometry'].centroid
            id_value = row[id]
            id_to_coords[id_value] = (centroid.y, centroid.x)  # (lat, lon)

        # Get the COMIDs from your self.forcing
        ds_ids = ds[id].values

        # Create coordinate arrays
        lats = np.array([id_to_coords[id_value][0] for id_value in ds_ids])
        lons = np.array([id_to_coords[id_value][1] for id_value in ds_ids])

        # Rename the COMID dimension to something more generic (we'll change it later)
        ds_renamed = ds.rename({id: 'location'})

        # Assign the new coordinates
        ds_renamed.coords['latitude'] = ('location', lats)
        ds_renamed.coords['longitude'] = ('location', lons)

        # Swap the dimension from COMID to lat/lon
        # First stack the lat/lon into a multi-index
        ds_stacked = ds_renamed.set_index(location=('latitude', 'longitude'))

        # Then unstack to create separate dimensions
        ds_final = ds_stacked.unstack('location')

        # resample the time dimension to daily frequency
        ds_final = ds_final.resample(time='1D').mean()

        # Assign the final self.forcing to self.forcing
        self.forcing = ds_final

        return

    def init_pet(self):
        """Initialize potential evapotranspiration (PET) calculation."""
        if self.pet_method != "hamon":
            raise ValueError(f"Invalid PET method: {self.pet_method}. Only 'hamon' is supported.")
        
        # FIXME: More flexibility is needed in the future for different PET methods
        # extract the latitude from the catchment object, first checking if `cat` has a crs
        if not self.cat.crs:
            warnings.warn(
                "Catchment (cat) does not have a CRS. A default of "
                "EPSG:4326 will be used for latitude calculation.",
                UserWarning)
            self.cat.set_crs(epsg=4326, inplace=True)
        # if crs is not set to EPSG:4326, then convert it to EPSG:4326
        elif self.cat.crs != 'EPSG:4326':
            self.cat.to_crs(epsg=4326, inplace=True)

        # calculate the latitude from the catchment geometry
        lats = self.forcing['latitude'].to_numpy()[0]

        # Calculate for each lat/lon point
        # Calculate PET while preserving dimensions
        temp_df = self.forcing['temp'].to_dataframe()
        temp_df = temp_df.reset_index().set_index('time')
        temp_df = temp_df.drop(columns=['latitude', 'longitude'], errors='ignore')

        pet = pyet.hamon(
            tmean=temp_df['temp'],
            lat=lats,
        )

        # assign pet to the self.forcing
        self.forcing['pet'] = pet

        # Set attributes for the new variable
        self.forcing['pet'].attrs = {
            'long_name': 'Potential Evapotranspiration',
            'units': 'millimeter / day',
            'source': 'Potential evapotranspiration estimated using Hamon, 1963, Trans. Am. Soc. Civ. Eng.'
        }

        # fix dimension orders
        # Create a spatial mask of ones with lat/lon dimensions
        spatial_ones = xr.ones_like(self.forcing['temp'].isel(time=0))

        # Multiply your time-series variable by the ones to broadcast it
        self.forcing['pet'] = self.forcing['pet'] * spatial_ones

        # Assuring the order of dimensions is time, latitude, longitude
        self.forcing['pet'] = self.forcing['pet'].transpose('time', 'latitude', 'longitude')

        return

    def init_streamflow(self):
        """Initialize streamflow data."""
        if self.streamflow is None:
            self.forcing['q_obs'] = xr.ones_like(self.forcing.temp)
        else:
            # adding streamflow data to the forcing object
            ds_obs = self.streamflow.sel(time=slice(self.settings['start_date'], self.settings['end_date']))
            # Create a spatial mask of ones with lat/lon dimensions
            spatial_ones = xr.ones_like(self.forcing.temp.isel(time=0))

            # Multiply your time-series variable by the ones to broadcast it
            self.forcing['q_obs'] = ds_obs * spatial_ones

        # adjust attributes for the streamflow data
        self.forcing['q_obs'].attrs = {
            'long_name': 'Mean observed daily discharge',
            'units': 'millimeter / day',
        }

        return

    def init_elev_bands(self):
        """Initialize elevation bands."""
        # if extra information is provided in the input files
        if self.elev_bands is not None:
            # read the elevation bands from the provided path
            if isinstance(self.elev_bands, (PathLike, str)):
                elev_bands_file = pd.read_csv(self.elev_bands, index_col=0, header=0)
                elev_bands_value = elev_bands_file.iloc[0, 0].values()[0]
            else:
                raise TypeError("elev_bands must be a PathLike or a string.")

        else:
            elev_bands_value = 1000

        # create a numpy.array of the elevation bands
        data = np.array([[[elev_bands_value]]])  # Shape: (elevation_band, latitude, longitude)

        # default prec_frac and area_frac values
        prec_frac = np.array([[[1.0]]])  # Shape: (elevation_band, latitude, longitude)
        area_frac = np.array([[[1.0]]])  # Shape: (elevation_band, latitude, longitude)

        ds = xr.Dataset(
            {
                'mean_elev': (['elevation_band', 'latitude', 'longitude'], data),
                'area_frac': (['elevation_band', 'latitude', 'longitude'], area_frac),
                'prec_frac': (['elevation_band', 'latitude', 'longitude'], prec_frac),
            },
            coords={
                'latitude': self.forcing.latitude, # Assuming latitude is defined in the forcing data
                'longitude': self.forcing.longitude, # Assuming longitude is defined in the forcing data
            }
        )

        # Assign the elevation_band as a coordinate variable
        ds = ds.assign_coords(elevation_band=np.arange(1, len(data) + 1))

        # Add attributes
        ds['mean_elev'].attrs = {
            'long_name': 'Mid-point elevation of each elevation band',
            'units': 'm asl'
        }
        ds['area_frac'].attrs = {
            'long_name': 'Fraction of the catchment covered by each elevation band',
            'units': 'dimensionless'
        }
        ds['prec_frac'].attrs = {
            'long_name': 'Fraction of catchment precipitation that falls on each elevation band - same as area_frac',
            'units': 'dimensionless'
        }
        ds['elevation_band'].attrs = {
            'long_name': 'elevation_band',
            'units': 'dimensionless'
        }
        ds['latitude'].attrs = {
            'long_name': 'latitude',
            'units': 'degreesN'
        }
        ds['longitude'].attrs = {
            'long_name': 'longitude',
            'units': 'degreesE'
        }
        ds['elevation_band'].attrs = {
            'long_name': 'elevation_band',
            'units': 'dimensionless'
        }

        self.elev_bands = ds

        return

    def init_model_file(
        self,
        base_path: str | PathLike, # type: ignore
        model_n: int,
    ) -> None:
        """Initialize the model file for the given model number."""

        paths_dict = {
            'settings': os.path.join(base_path, 'settings/'),
            'input': os.path.join(base_path, 'input/'),
            'output': os.path.join(base_path, 'output/'),
        }

        # if `start_evaluation` and `end_evaluation` are not provided in the settings,
        # then use the start and end dates from the forcing data
        if 'start_evaluation' not in self.settings or 'end_evaluation' not in self.settings:
            if 'time' not in self.forcing.dims:
                raise ValueError("Forcing data must have a 'time' dimension to extract start and end dates.")
            start_eval = self.forcing.time[-10]
            end_eval = self.forcing.time.max()[0]
        else:
            start_eval = self.settings['start_evaluation']
            end_eval = self.settings['end_evaluation']

        # start_date and end_date are mandatory settings
        start_date = self.settings['start_date']
        end_date = self.settings['end_date']

        date_dict = {
            'start_date': self._format_dates(start_date),
            'end_date': self._format_dates(end_date),
            'start_evaluation': self._format_dates(start_eval),
            'end_evaluation': self._format_dates(end_eval),
        }

        # create the content of the model file
        fm_catch = render_settings(
            paths=paths_dict,
            dates=date_dict,
            model=model_n
        )

        # return the rendered content
        return fm_catch
    
    def _format_dates(
        self,
        date: str,
    ) -> str:
        """Format the date string for the model file."""
        # read the date string and convert it to a datetime object
        date_obj = parser.parse(date)
        return date_obj.strftime("%Y-%m-%d")
