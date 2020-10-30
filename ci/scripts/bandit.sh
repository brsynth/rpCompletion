#!/bin/bash

source extras/.env


bandit -r ${PACKAGE}
