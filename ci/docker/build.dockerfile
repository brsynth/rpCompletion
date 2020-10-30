ARG IMAGE
FROM ${IMAGE}

RUN conda update -n base -c defaults conda

COPY recipe/conda_build_env.yaml /tmp/conda_build_env.yaml

RUN conda env create -n build --file /tmp/conda_build_env.yaml

ARG HOME
WORKDIR ${HOME}/ci

ENTRYPOINT ["conda", "run", "-n", "build"]
