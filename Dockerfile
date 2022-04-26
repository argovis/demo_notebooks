FROM conda/miniconda3

RUN apt-get update
RUN apt-get install -y nano gcc
RUN pip install notebook jupyter_contrib_nbextensions
RUN jupyter contrib nbextension install
RUN jupyter nbextension enable toc2/main
WORKDIR /books
COPY environment.yml /books/environment.yml
SHELL ["/bin/bash", "-c"]
RUN conda env create -f environment.yml && conda init bash && source ~/.bashrc && conda activate argovis_demos && conda install -y -c anaconda ipykernel conda-forge xarray dask netCDF4 bottleneck && python -m ipykernel install --user --name=argovis_demos && conda deactivate
CMD jupyter notebook --allow-root --ip=0.0.0.0
