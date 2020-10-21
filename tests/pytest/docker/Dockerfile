ARG IMAGE
FROM ${IMAGE}

RUN conda update -n base -c defaults conda

COPY tests/environment.yml        conda_env_run.yml
COPY tests/pytest/environment.yml conda_env_test.yml
ARG PKG
RUN conda env create -n test --file conda_env_run.yml \
 && conda env update -n test --file conda_env_test.yml

ARG HOME
WORKDIR ${HOME}/tests
ADD tests .

ADD ${PKG} ${PKG}

ENTRYPOINT ["conda", "run", "-n", "test"]
