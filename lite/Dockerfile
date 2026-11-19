# Build environment for the JupyterLite (runs-in-the-browser) version of the
# demo notebooks. This image only holds the build tools; the notebooks and the
# build script are mounted from the repository when the container runs.
#
# How to build the image, build the site and preview it: see lite/README.md.
#
# The versions pinned below decide which Pyodide (and so which Python, numpy,
# xarray, ...) the notebooks run on in the browser. jupyterlite-pyodide-kernel
# 0.8.7 uses Pyodide 314.0.6. Changing these pins is how the stack is upgraded.

FROM python:3.13-slim

RUN pip install --no-cache-dir \
    "jupyterlite-core==0.8.7" \
    "jupyterlite-pyodide-kernel[lock]==0.8.7" \
    "jupyter-server==2.21.1" \
    "nbconvert==7.17.2" \
    "build==1.6.1"

WORKDIR /work
CMD ["bash", "/repo/lite/build_site.sh"]
