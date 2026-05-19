# Argovis Demo Notebooks

This repository is a collection of Jupyter notebooks for learning how to search, download, and work with oceanographic data through [Argovis](https://argovis.colorado.edu/). The notebooks are organized into three folders:

## [`introduction/`](./introduction)

Start here. These notebooks introduce the core data types Argovis serves and the general patterns you'll use to query them:

- [`PROFILES.ipynb`](./introduction/PROFILES.ipynb) — working with point/profile datasets (Argo, GO-SHIP, etc.): querying by region and time, plotting maps and vertical profiles, and making T/S diagrams.
- [`GRIDS.ipynb`](./introduction/GRIDS.ipynb) — working with gridded products: querying a collection, subsetting in space, and loading results as an xarray.
- [`Argovis_JSON.ipynb`](./introduction/Argovis_JSON.ipynb) — for users working outside Python or in specialized applications, a tour of Argovis's raw JSON API responses.

## [`activities/`](./activities)

Hands-on, guided exercises that build on the introductory material. Each activity walks through a small scientific question end-to-end using the Argovis API.

- [`DOXY_ACTIVITY.ipynb`](./activities/DOXY_ACTIVITY.ipynb) — characterizing the vertical structure of dissolved oxygen in the ocean using BGC-Argo profiles and gridded products.

## [`dataset_specific_notebooks/`](./dataset_specific_notebooks)

Deep dives into individual datasets whose schemas or query patterns differ from the standard profile/grid examples covered in `introduction/`. Reach for these when you want to work with a particular product.

- [`Argo_trajectories.ipynb`](./dataset_specific_notebooks/Argo_trajectories.ipynb) — estimated Argo parking-depth trajectories and velocities.
- [`Argo_float_location_forecasts.ipynb`](./dataset_specific_notebooks/Argo_float_location_forecasts.ipynb) — the Argone float location forecast API.
- [`Intro_to_Atmospheric_Rivers.ipynb`](./dataset_specific_notebooks/Intro_to_Atmospheric_Rivers.ipynb) — atmospheric rivers, and Argovis's _extended objects_ schema for region-valued data.

---

## Running the notebooks

To run locally with Docker, clone the repo and mount it into the prebuilt image. First, install Git and Docker on your machine, then start the Docker application. Finally, run the following commands in the terminal:"

```
git clone https://github.com/argovis/demo_notebooks/
cd demo_notebooks
docker container run -p 8888:8888 -v $(pwd):/books argovis/notebooks jupyter notebook --allow-root --ip=0.0.0.0
```

On Windows, replace `$(pwd)` with the full path to your cloned `demo_notebooks` directory (e.g. `C:\Users\username\python_files\Argovis\demo_notebooks`) and run the command from PowerShell. In either case, copy the printed `http://127.0.0.1` URL into your browser to access the notebook environment.

You can refresh the image periodically with `docker image pull argovis/notebooks`, or update just the helper package from inside a notebook cell with `%pip install argovisHelpers`.

## See Also

- [Argovis' EarthCube 2022 submission](https://github.com/earthcube2022/ec22_mills_etal): illustrates dataset colocation, QC filtering, and interpolation.
