# cartopy normally downloads map data (coastlines etc.) from Natural Earth the
# first time a map is drawn. That download does not work in the browser, so
# the JupyterLite build (lite/build_site.sh) puts the files in the data/
# folder next to this file and ships them with the site.
#
# cartopy looks for a package with exactly this name, cartopy_userconfig, when
# it is imported, and calls update_config() from it. Nothing needs to import
# this package by hand.
from pathlib import Path


def update_config(config):
    # Tell cartopy to look in our data/ folder before trying to download.
    config['pre_existing_data_dir'] = Path(__file__).parent / 'data'
