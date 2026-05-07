# docker buildx build --platform linux/amd64,linux/arm64 -t argovis/notebooks --push .
FROM python:3.9

ENV PIP_CONFIG_FILE=/dev/null
RUN apt-get update
RUN apt-get install -y nano libproj-dev libgeos-dev libgdal-dev
RUN python -m pip install --upgrade pip "setuptools<81" wheel
RUN pip install jupyter numpy==1.22.3 requests==2.28.1 pandas==1.4.2 xarray==2022.11.0 matplotlib svgpath2mpl==1.0.0 scipy==1.9.3 geopandas==0.14.4
RUN pip install --no-cache-dir --no-build-isolation gsw
RUN pip install --no-cache-dir gsw-xarray
RUN pip install cartopy
RUN pip install argovisHelpers==0.0.35
RUN pip install netCDF4

WORKDIR /books
CMD jupyter notebook --allow-root --ip=0.0.0.0
