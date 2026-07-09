"""
Utilities for wedge-shaped Cartopy map plots with zebra boundaries.

Examples
--------
>>> fig, ax = wedge_map_boudary(lon_min=-25, lon_max=25)
>>> add_wedge_gridlines(ax)
>>> add_zebra_boundary(ax)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import Polygon
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.font_manager import FontProperties 

# ---------------------------------------------------------------------
# Core wedge map
# ---------------------------------------------------------------------

def wedge_map_boudary(ax: plt.Axes=None, projection: ccrs.Projection=None, lon_min: float=-25, lon_max: float=25, lat_min: float=75, lat_max: float=83, figsize: tuple =(7, 4), subplots=(1,1), land_color: str=None):
    """
    Create a wedge-shaped map with specified boundaries and optional land color.
    Parameters:
    - ax: Optional existing axis to draw on. If None, a new figure and axis will be created.
    - projection: Optional Cartopy projection. If None, AlbersEqualArea will be used with parameters based on the provided lat/lon limits.
    - lon_min, lon_max: Longitude limits of the map.
    - lat_min, lat_max: Latitude limits of the map.
    - figsize: Size of the figure if a new one is created.
    - land_color: Optional color for land features. If None, land will not be colored.
    Returns:
    - If ax is None, returns (fig, ax). Otherwise, returns ax.

    """

    # Create figure and axis if not provided
    if ax is None:
        # check if projection is provided, otherwise use AlbersEqualArea
        if projection is None:
            projection = ccrs.AlbersEqualArea(
                central_longitude=(lon_min+lon_max)/2,
                central_latitude=(lat_min+lat_max)/2,
                standard_parallels=(lat_min, lat_max)
            )
        fig, ax = plt.subplots(subplots[0], subplots[1], figsize=figsize, 
                                subplot_kw={'projection': projection}, constrained_layout=True)
    
    if isinstance(ax, np.ndarray):
        ax_list = ax.flat
    else:
        ax_list = [ax]

    for ax_i in ax_list:
        # Set the extent of the map
        ax_i.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())


        # Add land feature    
        if land_color is not None:
            ax_i.add_feature(cfeature.LAND, zorder=0, facecolor=land_color)
        
        # Draw boundary with all 4 edges properly defined
        n = 100
        bottom = [(lon, lat_min) for lon in np.linspace(lon_min, lon_max, n)]
        right  = [(lon_max, lat) for lat in np.linspace(lat_min, lat_max, n)]
        top    = [(lon, lat_max) for lon in np.linspace(lon_max, lon_min, n)]
        left   = [(lon_min, lat) for lat in np.linspace(lat_max, lat_min, n)]

        vertices = bottom + right + top + left
        boundary = Path(vertices)
        ax_i.set_boundary(boundary, transform=ccrs.PlateCarree())

    # return ax and fig if created, otherwise return ax
    if 'fig' in locals():
        return fig, ax
    else:
        return ax
    

# ---------------------------------------------------------------------
# Gridlines + labels
# ---------------------------------------------------------------------

def add_wedge_gridlines(ax, latitudes=[76, 77, 78, 79, 80, 81, 82],
                         longitudes=[-20, -10, 0, 10, 20],
                         lat_text=26, lon_text_bottom=74.8, lon_text_top=83.2,
                         zorder=11, fontsize=9,
                         show_top=True, show_bottom=True,
                         show_left=True, show_right=True,
                         lines=True,
                         font_path=None):
    
    """
    Add gridlines and custom labels to a wedge-shaped map.
    Parameters:
    - ax: The axis to draw on.
    - latitudes: List of latitudes for horizontal gridlines.
    - longitudes: List of longitudes for vertical gridlines.
    - lat_text: Distance from the edge for latitude labels.
    - lon_text_bottom: Latitude for longitude labels at the bottom edge.
    - lon_text_top: Latitude for longitude labels at the top edge.
    - zorder: Z-order for gridlines and labels.
    - fontsize: Font size for labels.
    - show_top, show_bottom, show_left, show_right: Booleans to control which edges show labels.
    - font_path: Path to a custom .otf or .ttf font file.
    """

    font_prop = FontProperties(fname=font_path) if font_path else None
    # Draw gridlines without labels
    if lines:
        gl = ax.gridlines(draw_labels=False, zorder=zorder, xlocs=longitudes, ylocs=latitudes)

    # write labels manually to control their position and formatting
    lat_sides = [(-lat_text, 'right', show_left), (lat_text, 'left', show_right)]
    lon_sides = [(lon_text_bottom, 'top', show_bottom), (lon_text_top, 'bottom', show_top)]  # adjust lat_max for top

    for x, ha, show in lat_sides:
        if show:
            for lat in latitudes:
                ax.text(x, lat, f'{lat}°N', transform=ccrs.PlateCarree(),
                        va='center', ha=ha, fontproperties=font_prop, fontsize=fontsize)

    for y, va, show in lon_sides:
        if show:
            for lon in longitudes:
                label = f'{abs(lon)}°{"W" if lon<0 else "E" if lon>0 else ""}'
                ax.text(lon, y, label, transform=ccrs.PlateCarree(),
                        va=va, ha='center', fontproperties=font_prop, fontsize=fontsize)
    
# ---------------------------------------------------------------------
# Zebra boundary
# ---------------------------------------------------------------------

def add_zebra_boundary(
    ax,
    lon_min=-25, lon_max=25,
    lat_min=75, lat_max=83,
    longitudes=None,
    latitudes=None,
    lw=8,
    outline_lw=1.5,
    zorder=10,
    n_interp=100,
    corners = True,
):
    """
    Draw a zebra-striped boundary around a cartopy map axes.

    The boundary is drawn along the four edges of the given lon/lat extent,
    with alternating black/white segments aligned to the provided graticule
    tick positions. Corners are filled with a solid grey patch to prevent gaps.

    Parameters
    ----------
    ax         : cartopy GeoAxes to draw on
    lon_min/max: longitude extent of the boundary
    lat_min/max: latitude extent of the boundary
    longitudes : interior longitude tick positions for stripe breaks
    latitudes  : interior latitude tick positions for stripe breaks
    lw         : line width (pts) of the zebra stripes
    outline_lw : extra width (pts) added to the black outline beneath each stripe
    zorder     : base zorder for the boundary lines
    n_interp   : number of interpolation points per segment (for curved projections)
    """

    if longitudes is None:
        longitudes = np.array([-20, -10, 0, 10, 20], dtype=float)
    if latitudes is None:
        latitudes = np.array([76, 77, 78, 79, 80, 81, 82], dtype=float)

    # Build break points by inserting interior ticks that fall within the extent
    lon_breaks = sorted({lon_min, *[l for l in longitudes if lon_min < l < lon_max], lon_max})
    lat_breaks = sorted({lat_min, *[l for l in latitudes  if lat_min < l < lat_max], lat_max})

    plate_crs = ccrs.PlateCarree()
    proj_crs  = ax.projection

    # -- Segment projection ---------------------------------------------------

    def project_segments(breaks, fixed_value, axis):
        """
        Project each consecutive pair of break points into axes coordinates.

        Returns a list of (xs, ys) arrays — one per segment — already
        transformed into the axes' native projection.
        """
        segments = []
        for start, end in zip(breaks, breaks[1:]):
            sample = np.linspace(start, end, n_interp)
            lons, lats = (sample, np.full_like(sample, fixed_value)) if axis == "lon" \
                    else (np.full_like(sample, fixed_value), sample)
            xy = proj_crs.transform_points(plate_crs, lons, lats)
            segments.append((xy[:, 0], xy[:, 1]))
        return segments

    # -- Edge drawing ---------------------------------------------------------

    def draw_striped_edge(breaks, fixed_value, axis):
        """
        Draw one edge of the zebra border as alternating black/white stripes,
        each sitting on top of a slightly thicker black outline.
        """
        segments = project_segments(breaks, fixed_value, axis)
        for i, (xs, ys) in enumerate(segments):
            stripe_color = "black" if i % 2 == 0 else "white"
            shared_kwargs = dict(
                transform=ax.transData,
                zorder=zorder,
                solid_capstyle="butt",
            )
            # Outline layer — drawn first, slightly wider
            ax.plot(xs, ys, color="black",    lw=lw + outline_lw, **shared_kwargs)
            # Stripe layer — drawn on top
            ax.plot(xs, ys, color=stripe_color, lw=lw, **(shared_kwargs | {"zorder": zorder + 1}))

    # Draw all four edges
    draw_striped_edge(lon_breaks, lat_min, "lon")  # bottom
    draw_striped_edge(lon_breaks, lat_max, "lon")  # top
    draw_striped_edge(lat_breaks, lon_min, "lat")  # left
    draw_striped_edge(lat_breaks, lon_max, "lat")  # right

    if corners:
        # -- Corners --------------------------------------------------------------
        # Cover the intersection of horizontal and vertical edges with a solid patch
        # sized slightly larger than the line width to avoid visible gaps.
        corner_size = 1.2 * _lw_to_data(ax, lw)
        corners = [
            (lon_min, lat_min),
            (lon_min, lat_max),
            (lon_max, lat_min),
            (lon_max, lat_max),
        ]
        for lon, lat in corners:
            cx, cy  = _projected_corner(ax, lon, lat)
            t_along_lon = _edge_tangent(ax, lon, lat, axis="lon")
            t_along_lat = _edge_tangent(ax, lon, lat, axis="lat")
            _draw_corner_patch(ax, cx, cy, t_along_lon, t_along_lat, corner_size)


# ---------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------

def _projected_corner(ax, lon, lat):
    """Return the projected (x, y) coordinate of a lon/lat corner point."""
    return ax.projection.transform_point(lon, lat, ccrs.PlateCarree())


def _edge_tangent(ax, lon, lat, dlon=0.01, dlat=0.01, axis="lon"):
    """
    Return a unit tangent vector at (lon, lat) along the given axis,
    in the axes' projected coordinate space.
    """
    plate = ccrs.PlateCarree()
    if axis == "lon":
        p0 = ax.projection.transform_point(lon,        lat,        plate)
        p1 = ax.projection.transform_point(lon + dlon, lat,        plate)
    else:
        p0 = ax.projection.transform_point(lon,        lat,        plate)
        p1 = ax.projection.transform_point(lon,        lat + dlat, plate)

    delta = np.array(p1) - np.array(p0)
    return delta / np.linalg.norm(delta)


def _lw_to_data(ax, lw_pts):
    """Convert a line width in points to an equivalent distance in data coordinates."""
    lw_px = lw_pts * ax.figure.dpi / 72
    inv   = ax.transData.inverted()
    x0    = inv.transform((0,     0))[0]
    x1    = inv.transform((lw_px, 0))[0]
    return abs(x1 - x0)


def _draw_corner_patch(ax, x, y, t1, t2, size, color="grey", zorder=200):
    """
    Draw a filled square patch at a boundary corner.

    The patch is oriented along the two edge tangent vectors (t1, t2) so it
    aligns correctly under any map projection, covering the join between the
    horizontal and vertical zebra edges.
    """
    half = size / 2
    v1, v2 = t1 * half, t2 * half

    # Four corners of the oriented square, offset from the centre point
    vertices = [
        (x - v1[0] - v2[0], y - v1[1] - v2[1]),
        (x + v1[0] - v2[0], y + v1[1] - v2[1]),
        (x + v1[0] + v2[0], y + v1[1] + v2[1]),
        (x - v1[0] + v2[0], y - v1[1] + v2[1]),
    ]

    ax.add_patch(
        Polygon(
            vertices,
            closed=True,
            facecolor=color,
            edgecolor="none",
            transform=ax.transData,
            zorder=zorder,
        )
    )