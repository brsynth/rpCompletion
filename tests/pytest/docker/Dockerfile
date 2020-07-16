ARG IMAGE
FROM ${IMAGE}

ARG PKG
COPY tests/pytest/requirements.txt requirements-test.txt
COPY ${PKG}/requirements.txt requirements-src.txt
# install requirements
RUN python3 -m pip install --upgrade pip \
 && python3 -m pip install --no-cache-dir --upgrade -r requirements-test.txt \
 && python3 -m pip install --no-cache-dir --upgrade -r requirements-src.txt

RUN conda install -y -c rdkit rdkit

ARG HOME
WORKDIR ${HOME}/tests

ADD tests .
ADD ${PKG} ${PKG}
