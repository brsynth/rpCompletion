ARG IMAGE
FROM ${IMAGE}

RUN conda update -n base -c defaults conda

COPY tests/pytest/conda_run.in  .
COPY tests/pytest/conda_test.in .
RUN conda install -y \
      --file conda_run.in \
      --file conda_test.in

ARG HOME
WORKDIR ${HOME}/tests
ADD tests .

ARG PKG
ADD ${PKG} ${PKG}
