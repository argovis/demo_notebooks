# cross compile for amd/arm:
# docker buildx create --name multiarch --driver docker-container --use
# docker buildx inspect --bootstrap
# docker buildx build --platform linux/amd64,linux/arm64 -t argovis/notebooks --push .

FROM python:3.14

ENV PIP_CONFIG_FILE=/dev/null
RUN apt-get update
RUN apt-get install -y nano libproj-dev libgeos-dev libgdal-dev
RUN python -m pip install --upgrade pip "setuptools<81" wheel
RUN pip install jupyter matplotlib svgpath2mpl==1.0.0 geopandas==0.14.4 cartopy netCDF4 argovisHelpers==0.0.41


WORKDIR /books
COPY pyproject.toml pyproject.toml
COPY notebook_helpers notebook_helpers
RUN pip install -e .
CMD jupyter notebook --allow-root --ip=0.0.0.0
