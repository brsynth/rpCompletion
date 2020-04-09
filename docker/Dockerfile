FROM brsynth/rpcache-advanced

COPY docker/requirements.txt .
# install requirements
RUN python -m pip install --upgrade pip \
 && python -m pip install -r requirements.txt

WORKDIR /home

COPY src .
