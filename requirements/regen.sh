#!/bin/sh

rm -f requirements/requirements.txt
cd docker
docker-compose --compatibility run --rm gen-req
