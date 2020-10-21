ARG IMAGE
FROM ${IMAGE}

RUN conda update -n base -c defaults conda

COPY tests/check/environment.yml conda_env_check.yml
ARG PKG
RUN conda env create -n check --file conda_env_check.yml

ARG HOME
WORKDIR ${HOME}
ADD tests tests
ADD ${PKG} ${PKG}

ENTRYPOINT ["conda", "run", "-n", "check"]
