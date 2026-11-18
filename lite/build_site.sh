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
"

cd /work
cp /repo/lite/jupyter_lite_config.json .

# 1. Package notebook_helpers as a wheel so the browser kernel can import it.
#    Build from a copy, because the repository is mounted read-only.
mkdir -p helpers_src
cp -r /repo/notebook_helpers /repo/pyproject.toml helpers_src/
python -m build --wheel --outdir dist helpers_src

# 2. Package the map data cartopy needs as a second wheel. cartopy cannot
#    download it in the browser, so it is downloaded here instead, from the
#    same address cartopy itself uses. See lite/cartopy_userconfig.
#    To add another Natural Earth dataset, repeat the download and unpack
#    lines with its address.
cp -r /repo/lite/cartopy_userconfig userconfig_src
MAP_DATA_DIR=userconfig_src/cartopy_userconfig/data/shapefiles/natural_earth/physical
mkdir -p "$MAP_DATA_DIR"
python -c "import urllib.request; urllib.request.urlretrieve('https://naturalearth.s3.amazonaws.com/50m_physical/ne_50m_coastline.zip', 'ne_50m_coastline.zip')"
python -m zipfile -e ne_50m_coastline.zip "$MAP_DATA_DIR"
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
