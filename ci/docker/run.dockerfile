ARG IMAGE
FROM ${IMAGE}

RUN conda update -n base -c defaults conda \
 && conda install -y -c conda-forge pyyaml

ARG HOME
WORKDIR ${HOME}/tests

ADD ci ci
ADD recipe recipe

RUN python ci/pytest/parse_recipe.py > /dev/null
RUN conda env create -n test --file ci/pytest/environment.yml

ARG PKG
ADD ${PKG} ${PKG}

ENV PACKAGE=$PKG
ENTRYPOINT ["conda", "run", "-n", "test", "python", "-m", "$PACKAGE"]
