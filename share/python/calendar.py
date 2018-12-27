#!/usr/bin/env python
# -*- coding: utf-8 -*-
# for '║'
import os
import sys
import subprocess
import numpy as np
import time
import datetime


ctime = time.localtime(time.time())
c_year  = ctime.tm_year
c_month = ctime.tm_mon
c_day   = ctime.tm_mday
c_hour  = ctime.tm_hour
c_min   = ctime.tm_min
c_sec   = ctime.tm_sec

s_datetime = datetime.datetime(2015,2,17,13,0,0)
c_datetime = datetime.datetime(c_year, c_month, c_day, c_hour, c_min, c_sec)

print (c_datetime - s_datetime).days
print (c_datetime - s_datetime).seconds
print (c_datetime - s_datetime).total_seconds()
