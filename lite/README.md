# JupyterLite build of the demo notebooks

This folder builds a version of the notebooks that runs entirely in the
browser, with no server and no install. The result is a folder of static files
that GitHub Pages serves.

Everything here runs in containers. All commands are run from the
**repository root**.

## Build and preview locally

1. Build the image that holds the build tools. Only needed once, and again if
   `lite/Dockerfile` changes.

   ```bash
   docker build -t argovis/lite-builder lite/
   ```

2. Build the site into `./_output` (ignored by git). Rerun this after any
   change to a notebook, a helper, or a file in this folder.

   ```bash
   docker run --rm -v $(pwd):/repo:ro -v $(pwd)/_output:/out argovis/lite-builder
   ```

3. Serve it, then open http://localhost:8000/ . The terminal should print
   `Serving on http://localhost:8000/ with browser caching switched off`.

   ```bash
   docker run --rm -p 8000:8000 -v $(pwd)/_output:/srv -v $(pwd)/lite/serve_site.py:/serve_site.py:ro -w /srv python:3.13-slim python /serve_site.py
   ```

The build needs internet access: it downloads Pyodide's package list, Python
packages, map data from Natural Earth, and two data files from GitHub.

## When the browser shows something stale

After a rebuild, reload the page and restart the kernel. If the old behaviour
is still there, it is one of these:

- **An old copy of a notebook.** JupyterLite saves every notebook you open
  into the browser's storage and shows that copy from then on, even after the
  site is rebuilt. In the JupyterLite file browser, right-click the notebook
  and choose Delete to get the site's version back.
- **Old cached files.** Use `lite/serve_site.py` (step 3 above), not
  `python -m http.server`, which lets browsers keep stale copies for up to an
  hour. To clear copies cached earlier: in Firefox, Cmd+Shift+Delete, tick only
  "Temporary cached files and pages".
- **When in doubt**, open the site in a private window. It starts with no
  cache and no saved notebooks.

## Memory

Every open notebook runs its own copy of Python in the browser, and closing a
notebook's tab does not stop it. After opening several notebooks, new ones
fail with out-of-memory or odd import errors. Use **Kernel → Shut Down All
Kernels**, or reload the page, and open one notebook at a time.

## Common changes

| To do this | Edit this |
|---|---|
| Publish another notebook | Add its path to `NOTEBOOKS` in `lite/build_site.sh` |
| Publish a data file beside a notebook | Commit the file, add its path to `DATA_FILES` in `lite/build_site.sh` |
| A helper in `notebook_helpers/` imports a new package | Add the package to `dependencies` in the root `pyproject.toml` |
| A notebook imports a package that is not part of Pyodide | Add it to `specs` in `lite/jupyter_lite_config.json` |
| A map needs another Natural Earth dataset | Add a `category/name` line to `MAP_DATASETS` in `lite/build_site.sh` |
| Upgrade Pyodide / Python / the scientific stack | Change the pinned versions in `lite/Dockerfile`, rebuild the image, retest every notebook |

## Rules for notebooks on this site

These come from how the browser kernel works. Breaking one gives a
"No module named ..." error.

- **Name a package in an `import` line, or list it in `pyproject.toml`.** The
  kernel loads a package when a cell imports it by name. A package that is
  only imported inside `notebook_helpers` is loaded because it is listed in
  the root `pyproject.toml`.
- **No `%magic` or `!shell` lines in a cell that has imports.** They stop the
  kernel from reading the cell's imports.
- **No `%pip install`.** Packages and their versions are fixed when the site
  is built.
- **Downloads must go through `requests`, to a server that allows requests
  from other websites.** `urllib` does not work in the browser. The Argovis
  API works; GitHub's raw file server does not, which is why the AR6 files are
  fetched at build time.

## What is in this folder

- `Dockerfile` — the build tools, with pinned versions. These pins decide
  which Pyodide the notebooks run on (currently Pyodide 314.0.6, Python 3.14).
- `build_site.sh` — the build steps, run inside the container. Lists the
  notebooks and data files to publish.
- `jupyter_lite_config.json` — packages to add to the site that are not part
  of Pyodide. `gsw` is listed by the full address of its WebAssembly file
  because the build tool does not find it by name.
- `cartopy_userconfig/` — a small package that carries the Natural Earth map
  data and tells cartopy to use it, because cartopy cannot download map data
  in the browser.
- `serve_site.py` — the local preview server.

## Publishing

`.github/workflows/lite.yml` runs the same image and build script on every
push to `main` and publishes the result to GitHub Pages. It can also be run by
hand from the repository's Actions tab.
