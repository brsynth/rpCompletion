#!/bin/bash

source extras/.env

# stop the build if there are Python syntax errors or undefined names
flake8 ${PACKAGE} --count --select=F82,E999 --show-source --statistics
