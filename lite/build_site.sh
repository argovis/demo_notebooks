#!/bin/bash
# Builds the JupyterLite site. Runs inside the container defined in
# lite/Dockerfile, which expects:
#   /repo  the demo_notebooks repository, mounted read-only
#   /out   an empty folder where the finished site is written
# All intermediate files stay in /work inside the container, so nothing in
# the repository is modified.
set -e

# Notebooks to publish, as paths relative to the repository root.
# To add a notebook to the site, add a line here.
NOTEBOOKS="
introduction/PROFILES.ipynb
introduction/Argovis_JSON.ipynb
introduction/GRIDS.ipynb
activities/ENSO_and_Equatorial_Kelvin_Waves.ipynb
activities/Mediterranean_Sea_Outflow.ipynb
activities/Meridional_Overturning_Circulation.ipynb
activities/Ocean_Stratification_Advanced.ipynb
activities/Oxygen_Vertical_Structure.ipynb
activities/Equatorial_Upwelling.ipynb
activities/Marine_Heatwaves.ipynb
activities/Ocean_Stratification_Intro.ipynb
activities/Subsurface_Ocean_Features_and_Water_Masses.ipynb
activities/Tropical_Cyclones_and_Ocean_Changes.ipynb
dataset_specific_notebooks/Argo_float_location_forecasts.ipynb
dataset_specific_notebooks/Intro_to_Atmospheric_Rivers.ipynb
"

cd /work
cp /repo/lite/jupyter_lite_config.json .

# 1. Package notebook_helpers as a wheel so the browser kernel can import it.
#    Build from a copy, because the repository is mounted read-only.
#    The dependency list in pyproject.toml is what makes the browser load
#    cartopy, jinja2 etc. along with notebook_helpers.
mkdir -p helpers_src
cp -r /repo/notebook_helpers /repo/pyproject.toml helpers_src/
python -m build --wheel --outdir dist helpers_src

# 2. Package the map data cartopy needs as a second wheel. cartopy cannot
#    download it in the browser, so it is downloaded here instead, from the
#    same address cartopy itself uses. See lite/cartopy_userconfig.
#
#    cartopy chooses a resolution from the size of the map being drawn
#    (110m for large areas, 50m under 50 degrees wide, 10m under 15 degrees),
#    so every dataset is bundled at all three resolutions.
#    To add another Natural Earth dataset, add a "category/name" line.
MAP_DATASETS="
physical/coastline
physical/land
cultural/admin_0_boundary_lines_land
"
MAP_RESOLUTIONS="110m 50m 10m"

cp -r /repo/lite/cartopy_userconfig userconfig_src
MAP_DATA_ROOT=userconfig_src/cartopy_userconfig/data/shapefiles/natural_earth
for dataset in $MAP_DATASETS; do
    category=$(dirname "$dataset")
    name=$(basename "$dataset")
    mkdir -p "$MAP_DATA_ROOT/$category"
    for resolution in $MAP_RESOLUTIONS; do
        zip_name="ne_${resolution}_${name}.zip"
        zip_url="https://naturalearth.s3.amazonaws.com/${resolution}_${category}/${zip_name}"
        python -c "import urllib.request; urllib.request.urlretrieve('$zip_url', '$zip_name')"
        python -m zipfile -e "$zip_name" "$MAP_DATA_ROOT/$category"
    done
done
python -m build --wheel --outdir dist userconfig_src

# 3. Copy the notebooks into the site's content folder, keeping their
#    subfolders, and clear their saved outputs to keep the downloads small.
for notebook in $NOTEBOOKS; do
    mkdir -p "content/$(dirname "$notebook")"
    cp "/repo/$notebook" "content/$notebook"
    jupyter nbconvert --clear-output --inplace "content/$notebook"
done

# 4. Build the static site. The packages named in jupyter_lite_config.json
#    are resolved here and recorded in the site's pyodide-lock.json.
#    Note: gsw is listed there by the full URL of its WebAssembly wheel on
#    PyPI. The resolver does not yet recognise that wheel's platform tag
#    (pyemscripten_2026_0) when gsw is requested by name alone.
#    The build deletes and recreates its output folder, which is not allowed
#    for a mounted folder like /out, so build inside the container first.
jupyter lite build --contents content --output-dir /work/_output

# 5. Replace whatever is in /out with the finished site.
find /out -mindepth 1 -delete
cp -r /work/_output/. /out/

echo "Site written to /out"
