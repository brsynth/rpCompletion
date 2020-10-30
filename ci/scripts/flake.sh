#!/bin/bash

source extras/.env


# stop the build if there are Python syntax errors or undefined names
flake8 ${PACKAGE} --count --select=E9,F63,F7,F82 --show-source --statistics
# exit-zero treats all errors as warnings. The GitHub editor is 127 chars wide
flake8 ${PACKAGE} --count --ignore=E272,E501,E266,E241,E226,E251,E303,E221 --exit-zero --max-complexity=10 --max-line-length=127 --statistics
# exit-zero treats all errors as warnings. The GitHub editor is 127 chars wide
flake8 tests --count --ignore=E272,E501,E266,E241,E226,E251,E303,E221,E122,E211,E302 --exit-zero --max-complexity=10 --max-line-length=127 --statistics
