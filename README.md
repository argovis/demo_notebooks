[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/argovis/demo_notebooks/HEAD)

# Argovis Demo Notebooks

This repository contains notebooks meant to be followed by individuals who want to learn more about searching, downloading and processing data with Argovis. Start with [Intro to 
Argovis](https://github.com/argovis/demo_notebooks/blob/main/Intro_to_Argovis.ipynb), and then follow along with any other notebooks that are relevant to your interests.

## Running the notebooks using Docker
Create a free Docker account at https://www.docker.com/ and install Docker on your laptop. Check that git is also installed on your laptop. With Docker installed (and running) on your laptop clone this repo, change directory into it, and mount it into a containerized environment:
```
git clone https://github.com/argovis/demo_notebooks/
cd demo_notebooks
docker container run -p 8888:8888 -v $(pwd):/books argovis/notebooks jupyter notebook --allow-root --ip=0.0.0.0
```
After a moment, several URLs will be printed to the terminal. Copy the one beginning with http://127.0.0.1 to your browser of choice to access the notebook environment. You may want to update your image from time to time with docker image pull argovis/notebooks, as we regularly update dependencies, especially https://pypi.org/project/argovisHelpers/⁠. Instead of updating the whole image, if you prefer you can just update the helper package, by running a notebook cell %pip install argovisHelpers.
![image](https://github.com/argovis/demo_notebooks/assets/14894641/6a115040-c70e-44cf-88b6-b12b9bf0f4c7)

## See Also

 - [Argovis' EarthCube 2022 submission](https://github.com/earthcube2022/ec22_mills_etal): illustrates dataset colocation, QC filtering, and interpolation.
