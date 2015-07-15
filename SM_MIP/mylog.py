#!/bin/env python
###################################
# Author: Jiang Li
# Email:  jiang_li@affymetrix.com
# Date:   Fri Nov 14 15:28:14 2014
###################################
import sys
import logging

logging.basicConfig(level=10,
        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
        datefmt='%a, %d %b %Y %H:%M:%S',
        stream=sys.stderr,
        filemode="w"
        )

def info(msg):
    logging.info(msg)

def debug(msg):
    logging.info(msg)

def warn(msg):
    logging.warning(msg)

def error(msg):
    logging.error(msg)



