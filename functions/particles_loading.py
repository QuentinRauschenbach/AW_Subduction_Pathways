import numpy as np
import xarray as xr
import pandas as pd
import datetime as dt
import gsw
from tqdm import tqdm


###
#
# Functions for loading and processing FESOM particle tracking output
#
###

##################################################################################################
#
# READ THE FILE
#
##################################################################################################

# Convert cftime to datetime
def noleap_to_datetime(cftime_dt):
    return dt.datetime(
        cftime_dt.year,
        cftime_dt.month,
        cftime_dt.day,
        cftime_dt.hour,
        cftime_dt.minute,
        cftime_dt.second,
    )

# Read File
def read_and_convert_time(file):
    """
    Reads a NetCDF file using xarray and attempts to decode the time variable. If decoding fails, it falls back to reading the time as raw values and converting them to datetime objects.

    Parameters:
    file (str): The path to the NetCDF file to be read.

    Returns:
    xarray.Dataset: The dataset with the time variable converted to datetime objects.
    """
    
    try:
        ds = xr.open_dataset(file, decode_times=True, )

        ds = ds.assign_coords(time=pd.DatetimeIndex([noleap_to_datetime(t) for t in ds.time.values]))
    except:
        print("Failed to decode times, opening without decoding.")
        ds = xr.open_dataset(file, decode_times=False)

        seconds = ds.time.values.astype("float64")

        seconds[seconds > 1e30] = np.nan

        origin = pd.Timestamp("2006-01-01")
        time = origin + pd.to_timedelta(seconds, unit="s")

        ds = ds.assign_coords(time=time)
    return ds

# Read Meta data from global attrs
def extract_metadata(ds, print_metadata=True):
    """
    Extracts metadata from the dataset attributes and calculates the release date.

    Parameters:
    ds (xarray.Dataset): The dataset from which to extract metadata.

    Returns:
    """

    STEPS_PER_DAY  = int(ds.attrs["steps_per_day"] / ds.attrs["n_out"]) # output steps per day
    DOY            = ds.attrs["release_start_doy"].item()
    YEAR           = ds.attrs["release_start_year"].item()

    # calculate release date
    release_date = pd.Timestamp(year=int(YEAR), month=1, day=1) + pd.Timedelta(days=int(DOY) - 1)

    print(f"Release date: {release_date}")
    print(f" - Day of year: {DOY}")
    print(f" - Out steps per day: {STEPS_PER_DAY}")

    return STEPS_PER_DAY, DOY, YEAR, release_date

##################################################################################################
#
# CALCULATE BASIC METRICS AND ADD TO THE DATASET
#
##################################################################################################

def calc_vorticity_strain_ow(du_dx, du_dy, dv_dx, dv_dy):
    """
    Calculate relative vorticity, strain rate, and Okubo-Weiss parameter
    from horizontal velocity gradients.

    Parameters
    ----------
    du_dx : array-like
        Partial derivative of u with respect to x (s^-1)
    du_dy : array-like
        Partial derivative of u with respect to y (s^-1)
    dv_dx : array-like
        Partial derivative of v with respect to x (s^-1)
    dv_dy : array-like
        Partial derivative of v with respect to y (s^-1)

    Returns
    -------
    vorticity : np.ndarray
        Relative vorticity ζ = dv/dx - du/dy (s^-1)
        Positive = cyclonic (Northern Hemisphere)
    strain : np.ndarray
        Total strain rate S = sqrt(S_n^2 + S_s^2) (s^-1)
        S_n = du/dx - dv/dy  (normal strain)
        S_s = dv/dx + du/dy  (shear strain)
    okubo_weiss : np.ndarray
        Okubo-Weiss parameter W = S^2 - ζ^2 (s^-2)
        W > 0 : strain-dominated (fronts, filaments)
        W < 0 : vorticity-dominated (eddy cores)
        W ~ 0 : background flow
    """
    du_dx = np.asarray(du_dx)
    du_dy = np.asarray(du_dy)
    dv_dx = np.asarray(dv_dx)
    dv_dy = np.asarray(dv_dy)

    # Relative vorticity
    vorticity = dv_dx - du_dy

    # Strain components
    S_normal = du_dx - dv_dy   # normal strain (stretching)
    S_shear  = dv_dx + du_dy   # shear strain

    # Total strain rate magnitude
    strain = np.sqrt(S_normal**2 + S_shear**2)

    # Okubo-Weiss parameter
    okubo_weiss = S_normal**2 + S_shear**2 - vorticity**2

    # Horizontal Divergence
    # Negative values = Convergence (downward forcing)
    # Positive values = Divergence (upwelling forcing)
    divergence = du_dx + dv_dy

    return vorticity, strain, okubo_weiss, divergence


def calc_rossby(vorticity, lat):
    """
    Calculate Rossby number from relative vorticity and latitude.

    Parameters
    ----------
    vorticity : array-like
        Relative vorticity ζ (s^-1)
    lat : array-like
        Latitude in degrees

    Returns
    -------
    Ro : np.ndarray
        Rossby number Ro = |ζ| / f
    f : np.ndarray
        Coriolis parameter f = 2Ω sin(lat) (s^-1)
    """
    OMEGA = 2 * np.pi / 86164  # Earth rotation rate (rad/s)

    lat = np.asarray(lat)
    vorticity = np.asarray(vorticity)

    f = 2 * OMEGA * np.sin(np.deg2rad(lat))
    Ro = vorticity / np.abs(f)

    return Ro, f

def calc_density(temp, salt, depth, lat):
    """
    Calculate in-situ density and sigma0 from temperature, salinity, depth, and latitude using TEOS-10 (gsw).
    Parameters
    ----------
    temp : array-like
        potential temperature (°C)
    salt : array-like
        Practical salinity    depth : array-like
        Depth in metres (positive downward)
    lat : array-like
        Latitude in degrees

    returns
    -------
    rho : np.ndarray
        In-situ density (kg m^-3)
    sigma0 : np.ndarray
        Potential density anomaly referenced to the surface (kg m^-3)
    """
    pressure = gsw.p_from_z(-depth, lat)
    SA = gsw.SA_from_SP(salt, pressure, np.zeros_like(temp), lat)
    CT = gsw.CT_from_pt(SA, temp)
    rho = gsw.rho(SA, CT, pressure)
    sigma0 = gsw.sigma0(SA, CT)
    
    return rho, sigma0

def calc_density_column_old_no_chunks(temp_column, salt_column, depth_levels, lat):
    temp_column    = np.asarray(temp_column)    # (depth, drifter, time)
    salt_column    = np.asarray(salt_column)
    depth_levels   = np.asarray(depth_levels)   # (depth,)
    lat            = np.asarray(lat)            # (drifter, time)

    # Convert depth → pressure (dbar) for gsw
    pressure_levels = gsw.p_from_z(-depth_levels, 79.0)   # (nz,)

    # Unpack dimensions
    nz, nd, nt = temp_column.shape
    print(f"nz={nz}, nd={nd}, nt={nt}")

    # Transpose to (drifter, time, depth) then flatten to (drifter*time, depth)
    T_flat   = temp_column.transpose(1, 2, 0).reshape(nd*nt, nz)
    S_flat   = salt_column.transpose(1, 2, 0).reshape(nd*nt, nz)
    lat_flat = lat.reshape(nd*nt, 1) * np.ones((1, nz))
    lon_flat = np.zeros_like(lat_flat)

    print("T_flat shape:", T_flat.shape)          # should be (nd*nt, nz)
    print("T_flat range:", np.nanmin(T_flat), np.nanmax(T_flat))

    # valid mask - exclude terminated particles
    valid = ~np.isnan(lat.reshape(nd*nt))

    print(f"valid points: {valid.sum()} / {len(valid)}")

    # only compute gsw where valid
    SA  = np.full((nd*nt, nz), np.nan)
    CT  = np.full((nd*nt, nz), np.nan)
    rho = np.full((nd*nt, nz), np.nan)
    sigma0 = np.full((nd*nt, nz), np.nan)

    SA[valid]  = gsw.SA_from_SP(S_flat[valid], pressure_levels[np.newaxis, :],
                                lon_flat[valid], lat_flat[valid])
    CT[valid]  = gsw.CT_from_pt(SA[valid], T_flat[valid])
    rho[valid] = gsw.rho(SA[valid], CT[valid], pressure_levels[np.newaxis, :])
    sigma0[valid] = gsw.sigma0(SA[valid], CT[valid])

    # Reshape back to (depth, drifter, time)
    rho = rho.reshape(nd, nt, nz).transpose(2, 0, 1)
    sigma0 = sigma0.reshape(nd, nt, nz).transpose(2, 0, 1)

    # set all points to NaN where the temperature is larger than 30 
    invalid_temp = temp_column > 30
    rho[invalid_temp] = np.nan
    sigma0[invalid_temp] = np.nan
    print("rho shape:", rho.shape)                # should be (nz, nd, nt)
    print("sample rho profile:", rho[:5, 0, 0])  # should increase ~1026→1028
    #print("fraction negative drho_dz will tell us sign convention is correct")

    return rho, sigma0

### TEST ####
def calc_mld(sigma0_chunk, depth_levels, threshold=0.03):
    """
    sigma0_chunk : (nz, nc, nt) — already computed for this chunk
    depth_levels : (nz,)        — positive downward, metres
    returns mld  : (nc, nt)     — in metres, NaN where not found
    """
    nz, nc, nt = sigma0_chunk.shape
    mld = np.full((nc, nt), np.nan, dtype="float32")

    # Surface density = shallowest non-NaN level
    surf = sigma0_chunk[0, :, :]   # (nc, nt)

    for k in range(1, nz):
        exceeded   = (sigma0_chunk[k, :, :] - surf) > threshold  # (nc, nt)
        not_filled = np.isnan(mld)
        mld = np.where(exceeded & not_filled, depth_levels[k], mld)

    return mld

def calc_density_column(ds, chunk_size=200, mld_threshold=0.03):
    pres         = gsw.p_from_z(-np.asarray(ds.depth_levels), 79.0)
    depth_levels = np.asarray(ds.depth_levels)
    nz           = len(pres)
    nd           = ds.temp_column.shape[1]
    nt           = ds.temp_column.shape[2]

    ds["density_column"] = xr.DataArray(
        np.full((nz, nd, nt), np.nan, dtype="float32"),
        coords=ds["temp_column"].coords,
        dims=ds["temp_column"].dims,
        attrs={"long_name": "In-situ density column", "units": "kg m-3"}
    )
    ds["sigma0_column"] = xr.DataArray(
        np.full((nz, nd, nt), np.nan, dtype="float32"),
        coords=ds["temp_column"].coords,
        dims=ds["temp_column"].dims,
        attrs={"long_name": "Potential density anomaly column", "units": "kg m-3"}
    )
    # MLD is 2D — no depth dimension
    ds["mixed"] = xr.DataArray(
        np.full((nd, nt), np.nan, dtype="float32"),
        coords={k: ds["temp_column"].coords[k] 
                for k in ds["temp_column"].dims[1:]},  # drifter, time
        dims=ds["temp_column"].dims[1:],
        attrs={"long_name": "Mixed layer depth", "units": "m",
               "mld_threshold": f"delta_sigma0 > {mld_threshold} kg/m3"}
    )

    for start in tqdm(range(0, nd, chunk_size), desc="...Calculating column density in chunks"):
        end = min(start + chunk_size, nd)
        nc  = end - start

        T_chunk   = ds.temp_column[:, start:end, :].values.transpose(1, 2, 0).reshape(nc * nt, nz)
        S_chunk   = ds.salt_column[:, start:end, :].values.transpose(1, 2, 0).reshape(nc * nt, nz)
        lat_chunk = ds.lat[start:end, :].values.reshape(nc * nt, 1) * np.ones((1, nz))
        lon_chunk = np.zeros_like(lat_chunk)

        valid = ~np.isnan(ds.lat[start:end, :].values.reshape(nc * nt))

        SA_c     = np.full((nc * nt, nz), np.nan, dtype="float32")
        CT_c     = np.full((nc * nt, nz), np.nan, dtype="float32")
        rho_c    = np.full((nc * nt, nz), np.nan, dtype="float32")
        sigma0_c = np.full((nc * nt, nz), np.nan, dtype="float32")

        if valid.any():
            SA_c[valid]     = gsw.SA_from_SP(S_chunk[valid], pres[np.newaxis, :],
                                              lon_chunk[valid], lat_chunk[valid])
            CT_c[valid]     = gsw.CT_from_pt(SA_c[valid], T_chunk[valid])
            rho_c[valid]    = gsw.rho(SA_c[valid], CT_c[valid], pres[np.newaxis, :])
            sigma0_c[valid] = gsw.sigma0(SA_c[valid], CT_c[valid])

        invalid_temp = T_chunk.reshape(nc, nt, nz).transpose(2, 0, 1) > 30
        rho_c_3d     = rho_c.reshape(nc, nt, nz).transpose(2, 0, 1)
        sigma0_c_3d  = sigma0_c.reshape(nc, nt, nz).transpose(2, 0, 1)
        rho_c_3d[invalid_temp]    = np.nan
        sigma0_c_3d[invalid_temp] = np.nan

        ds["density_column"].values[:, start:end, :]  = rho_c_3d
        ds["sigma0_column"].values[:, start:end, :]   = sigma0_c_3d
        ds["mixed"].values[start:end, :]                = calc_mld(
            sigma0_c_3d, depth_levels, threshold=mld_threshold
        )

        #print(f"  chunk {start}:{end} done")

    return ds


def add_diagnostics_to_dataset(ds, lat_var="lat", calc_column=False):
    """
    Calculate submesoscale diagnostics and add them to an xarray Dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset containing du_dx, du_dy, dv_dx, dv_dy as variables
        and a latitude variable (default: 'lat')
    lat_var : str
        Name of the latitude variable in ds

    Returns
    -------
    ds : xarray.Dataset
        Original dataset with added variables:
        - density     : in-situ density ρ (kg m^-3)
        - sigma0      : potential density anomaly σ0 (kg m^-3)
        - vorticity   : relative vorticity ζ (s^-1)
        - strain      : total strain rate S (s^-1)
        - okubo_weiss : Okubo-Weiss parameter W (s^-2)
        - rossby      : Rossby number Ro (-)
        - f           : Coriolis parameter (s^-1)
    """

    density, sigma0 = calc_density(ds.temp, ds.salt, ds.depth, ds.lat)
    if calc_column:
        if "temp_column" in ds and "salt_column" in ds and "depth_levels" in ds:
            #density_column, sigma0_column = calc_density_column(ds.temp_column, ds.salt_column, ds.depth_levels, ds.lat, chunk_size=500)
            ds = calc_density_column(ds, chunk_size=200)

    vorticity, strain, okubo_weiss, divergence = calc_vorticity_strain_ow(
        ds["du_dx"], ds["du_dy"], ds["dv_dx"], ds["dv_dy"]
    )
    rossby, f = calc_rossby(vorticity, ds[lat_var])

    ds["density"] = xr.DataArray(
        density,
        coords=ds["temp"].coords,
        dims=ds["temp"].dims,
        attrs={"long_name": "In-situ density", "units": "kg m-3"}
    )
    ds["sigma0"] = xr.DataArray(
        sigma0,
        coords=ds["temp"].coords,
        dims=ds["temp"].dims,
        attrs={"long_name": "Potential density anomaly referenced to surface", "units": "kg m-3"}
    )
    if False:
        if "temp_column" in ds and "salt_column" in ds and "depth_levels" in ds:
            ds["density_column"] = xr.DataArray(
                density_column,
                coords=ds["temp_column"].coords,
                dims=ds["temp_column"].dims,
                attrs={"long_name": "In-situ density column", "units": "kg m-3"}
            )
            ds["sigma0_column"] = xr.DataArray(
                sigma0_column,
                coords=ds["temp_column"].coords,
                dims=ds["temp_column"].dims,
                attrs={"long_name": "Potential density anomaly column", "units": "kg m-3"}
            )


    ds["vort"] = xr.DataArray(
        vorticity,
        coords=ds["du_dx"].coords,
        dims=ds["du_dx"].dims,
        attrs={"long_name": "Relative vorticity", "units": "s-1"}
    )
    ds["strain"] = xr.DataArray(
        strain,
        coords=ds["du_dx"].coords,
        dims=ds["du_dx"].dims,
        attrs={"long_name": "Total strain rate", "units": "s-1"}
    )
    ds["divergence"] = xr.DataArray(
        divergence,
        coords=ds["du_dx"].coords,
        dims=ds["du_dx"].dims,
        attrs={"long_name": "Horizontal divergence", "units": "s-1"}
    )
    ds["okubo_weiss"] = xr.DataArray(
        okubo_weiss,
        coords=ds["du_dx"].coords,
        dims=ds["du_dx"].dims,
        attrs={"long_name": "Okubo-Weiss parameter", "units": "s-2"}
    )
    ds["Ro"] = xr.DataArray(
        rossby,
        coords=ds["du_dx"].coords,
        dims=ds["du_dx"].dims,
        attrs={"long_name": "Rossby number", "units": "-"}
    )
    ds["f"] = xr.DataArray(
        f,
        coords=ds[lat_var].coords,
        dims=ds[lat_var].dims,
        attrs={"long_name": "Coriolis parameter", "units": "s-1"}
    )

    return ds


def calc_drho_dz_and_Ri(temp_column, salt_column, depth_levels, particle_depth,
                        du_dz, dv_dz, lat):
    """
    Calculate vertical density gradient and Richardson number at the particle
    depth, using the layer above and below via centred differencing.

    Uses TEOS-10 (gsw) for density. Depth is used as pressure proxy under
    the Boussinesq / incompressible assumption used in FESOM.

    Parameters
    ----------
    temp_column : array-like, shape (n_drifter, n_time, n_depth)
        In-situ temperature (°C) of the water column at particle position
    salt_column : array-like, shape (n_drifter, n_time, n_depth)
        Practical salinity of the water column at particle position
    depth_levels : array-like, shape (n_depth,)
        Depth of each level in metres (positive downward)
    particle_depth : array-like, shape (n_drifter, n_time)
        Current depth of the particle in metres (positive downward)
    du_dz : array-like, shape (n_drifter, n_time)
        Vertical shear of u at particle position (s^-1)
    dv_dz : array-like, shape (n_drifter, n_time)
        Vertical shear of v at particle position (s^-1)
    lat : array-like, shape (n_drifter,) or (n_drifter, n_time)
        Latitude in degrees

    Returns
    -------
    drho_dz : np.ndarray, shape (n_drifter, n_time)
        Vertical density gradient dρ/dz (kg m^-4), negative = stable
    N2 : np.ndarray, shape (n_drifter, n_time)
        Buoyancy frequency squared N² = -(g/ρ₀) dρ/dz (s^-2)
    Ri : np.ndarray, shape (n_drifter, n_time)
        Richardson number Ri = N² / (du_dz² + dv_dz²)
        Ri < 0.25 : shear instability likely
        Ri < 1    : submesoscale-relevant
    """
    G    = 9.81       # m s^-2
    RHO0 = 1025.0     # reference density kg m^-3

    temp_column    = np.asarray(temp_column)    # (depth, drifter, time)
    salt_column    = np.asarray(salt_column)
    depth_levels   = np.asarray(depth_levels)   # (depth,)
    particle_depth = np.asarray(particle_depth)           # (drifter, time)
    du_dz          = np.asarray(du_dz)           # (drifter, time)
    dv_dz          = np.asarray(dv_dz)           # (drifter, time)
    lat            = np.asarray(lat)             # (drifter, time)

    # Convert depth → pressure (dbar) for gsw
    pressure_levels = gsw.p_from_z(-depth_levels, 79.0)   # (nz,)

    # Unpack dimensions
    nz, nd, nt = temp_column.shape
    print(f"nz={nz}, nd={nd}, nt={nt}")

    # Transpose to (drifter, time, depth) then flatten to (drifter*time, depth)
    T_flat   = temp_column.transpose(1, 2, 0).reshape(nd*nt, nz)
    S_flat   = salt_column.transpose(1, 2, 0).reshape(nd*nt, nz)
    lat_flat = lat.reshape(nd*nt, 1) * np.ones((1, nz))
    lon_flat = np.zeros_like(lat_flat)

    print("T_flat shape:", T_flat.shape)          # should be (nd*nt, nz)
    print("T_flat range:", np.nanmin(T_flat), np.nanmax(T_flat))

    # valid mask - exclude terminated particles
    valid = ~np.isnan(lat.reshape(nd*nt))
    print(f"valid points: {valid.sum()} / {len(valid)}")

    # only compute gsw where valid
    SA  = np.full((nd*nt, nz), np.nan)
    CT  = np.full((nd*nt, nz), np.nan)
    rho = np.full((nd*nt, nz), np.nan)

    SA[valid]  = gsw.SA_from_SP(S_flat[valid], pressure_levels[np.newaxis, :],
                                lon_flat[valid], lat_flat[valid])
    CT[valid]  = gsw.CT_from_t(SA[valid], T_flat[valid], pressure_levels[np.newaxis, :])
    rho[valid] = gsw.rho(SA[valid], CT[valid], pressure_levels[np.newaxis, :])

    # Reshape back to (depth, drifter, time)
    rho = rho.reshape(nd, nt, nz).transpose(2, 0, 1)
    print("rho shape:", rho.shape)                # should be (nz, nd, nt)
    print("sample rho profile:", rho[:5, 0, 0])  # should increase ~1026→1028
    print("fraction negative drho_dz will tell us sign convention is correct")

    # Find the depth level index closest to each particle
    # depth_levels: (nz,), particle_depth: (nd, nt)
    depth_diff = np.abs(depth_levels[:, np.newaxis, np.newaxis] -
                        particle_depth[np.newaxis, :, :])  # (nz, nd, nt)
    k = np.argmin(depth_diff, axis=0)  # (nd, nt)

    # Clip so we always have a level above and below
    k = np.clip(k, 1, nz - 2)

    # Centred difference: dρ/dz using layers above (k-1) and below (k+1)
    nd_idx = np.arange(nd)[:, np.newaxis]
    nt_idx = np.arange(nt)[np.newaxis, :]

    rho_above = rho[k - 1, nd_idx, nt_idx]
    rho_below = rho[k + 1, nd_idx, nt_idx]
    dz        = depth_levels[k + 1] - depth_levels[k - 1]  # (nd, nt)

    # dρ/dz positive downward
    drho_dz = (rho_below - rho_above) / dz  # (nd, nt)

    # N² = (g/ρ₀) * dρ/dz  (positive = stable)
    N2 = (G / RHO0) * drho_dz

    print("N2 range:", np.nanmin(N2), np.nanmax(N2))
    print("fraction negative N2:", np.nanmean(N2 < 0))  # should be <5%

    # Richardson number
    shear2 = du_dz**2 + dv_dz**2
    with np.errstate(divide="ignore", invalid="ignore"):
        Ri = np.where(shear2 > 1e-10, N2 / shear2, np.nan)

    print("Ri shape:", Ri.shape)
    print("Ri range:", np.nanmin(Ri), np.nanmax(Ri))
    print("fraction Ri < 0.25:", np.nanmean(Ri < 0.25))
    print("fraction Ri < 1:", np.nanmean(Ri < 1.0))

    return drho_dz, N2, Ri


def add_stratification_to_dataset(ds,
                                  depth_var="depth",
                                  lat_var="lat",
                                  du_dz_var="du_dz",
                                  dv_dz_var="dv_dz"):
    """
    Calculate drho_dz, N², and Ri and add them to an xarray Dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        Must contain: temp_column, salt_column, du_dz, dv_dz,
        and the particle depth and latitude variables.
    depth_levels : array-like
        Depth of each level in metres (positive downward)
    depth_var : str
        Variable name for particle depth in ds (default: 'depth')
    lat_var : str
        Variable name for latitude in ds (default: 'lat')
    du_dz_var, dv_dz_var : str
        Variable names for vertical velocity shear in ds

    Returns
    -------
    ds : xarray.Dataset
        Dataset with added drho_dz, N2, and Ri variables
    """
    import xarray as xr

    drho_dz, N2, Ri = calc_drho_dz_and_Ri(
        temp_column    = ds["temp_column"].values,
        salt_column    = ds["salt_column"].values,
        depth_levels   = ds.depth_levels.values,
        particle_depth = ds[depth_var].values,
        du_dz          = ds[du_dz_var].values,
        dv_dz          = ds[dv_dz_var].values,
        lat            = ds[lat_var].values,
    )

    # Use (drifter, time) dims/coords from an existing 2D variable
    ref     = ds[depth_var]
    coords  = ref.coords
    dims    = ref.dims

    ds["drho_dz"] = xr.DataArray(
        drho_dz, coords=coords, dims=dims,
        attrs={"long_name": "Vertical density gradient", "units": "kg m-4"}
    )
    ds["N2"] = xr.DataArray(
        N2, coords=coords, dims=dims,
        attrs={"long_name": "Buoyancy frequency squared", "units": "s-2"}
    )
    ds["Ri"] = xr.DataArray(
        Ri, coords=coords, dims=dims,
        attrs={"long_name": "Richardson number", "units": "-"}
    )

    return ds


##################################################################################################
#
# COMBINE
#
##################################################################################################

def read_particle_data(file, calc_column=False):
    """
    Reads particle data from a NetCDF file and extracts metadata.

    Parameters:
    file (str): The path to the NetCDF file to be read.

    Returns:
    xarray.Dataset: The dataset containing particle data with time converted to datetime objects.
    int: Steps per day extracted from the dataset attributes.
    int: Day of year extracted from the dataset attributes.
    int: Year extracted from the dataset attributes.
    pd.Timestamp: Release date calculated from the day of year and year.
    """

    ds = read_and_convert_time(file)
    
    try:
        STEPS_PER_DAY, DOY, YEAR, release_date = extract_metadata(ds)
    except KeyError as e:
        print(f"Metadata key error: {e}. Check dataset attributes.")
        STEPS_PER_DAY, DOY, YEAR, release_date = None, None, None, None
    
    try:
        ds = add_diagnostics_to_dataset(ds, calc_column=calc_column)
    except KeyError as e:
        print(f"Diagnostic calculation error: {e}. Check required variables for diagnostics.")
    return ds, STEPS_PER_DAY, DOY, YEAR, release_date
