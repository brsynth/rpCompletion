ARG IMAGE
FROM ${IMAGE}

RUN conda update -n base -c defaults conda

COPY ci/check/environment.yml conda_env_check.yml
RUN conda env create -n check --file conda_env_check.yml

ARG HOME
WORKDIR ${HOME}

ENTRYPOINT ["conda", "run", "-n", "check"]
