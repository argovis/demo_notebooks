[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/argovis/demo_notebooks/HEAD)

# Argovis Demo Notebooks

This repository contains notebooks meant to be followed by individuals who want to learn more about searching, downloading and processing data with Argovis. Start with [Intro to 
Argovis](https://github.com/argovis/demo_notebooks/blob/main/Intro_to_Argovis.ipynb), and then follow along with any other notebooks that are relevant to your interests.

## Running the notebooks using Docker
Create a free Docker account at https://www.docker.com/ and install Docker (https://www.docker.com/products/docker-desktop/) on your laptop. Check that git is also installed on your laptop. Before moving forward, please make sure no other jupyter notebooks are running on your laptop. With Docker installed (and running) on your laptop clone the `https://github.com/argovis/demo_notebooks/` repo, change directory into it, and mount it into a containerized environment:

```
git clone https://github.com/argovis/demo_notebooks/
cd demo_notebooks
docker container run -p 8888:8888 -v $(pwd):/books argovis/notebooks jupyter notebook --allow-root --ip=0.0.0.0
```

After a moment, several URLs will be printed to the terminal. Copy the one beginning with http://127.0.0.1 to your browser of choice to access the notebook environment. You may want to update your image from time to time with ```docker image pull argovis/notebooks```, as we regularly update dependencies, especially [argovisHelpers](https://pypi.org/project/argovisHelpers/). Instead of updating the whole image, if you prefer you can just update the helper package, by running a notebook cell ```%pip install argovisHelpers```.

## Edits to the above, if you are using Windows
Once you have cloned the demo_notebooks repository:
- use `Windows powershell` to `cd` to the `demo_notebooks` directory
- use `Windows powershell` to run the this command:

  ```
  docker container run -p 8888:8888 -v path_to_files:/books argovis/notebooks jupyter notebook --allow-root --ip=0.0.0.0
  ```

  Please note that `path_to_files` is the full path to the `demo_notebooks` directory on your machine, e.g. `C:\Users\username\python_files\Argovis\demo_notebooks`.
- After a moment, several URLs will be printed to the terminal. Copy the one beginning with http://127.0.0.1 to your browser of choice to access the notebook environment. 

## See Also

 - [Argovis' EarthCube 2022 submission](https://github.com/earthcube2022/ec22_mills_etal): illustrates dataset colocation, QC filtering, and interpolation.
