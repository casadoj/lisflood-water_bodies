#!/usr/bin/python3
# -*- coding: utf-8 -*-
""" Version: 1.11 """
""" Author: Christian Schwatke <christian.schwatke@tum.de>"""
""" Last change: 2022-01-10 """

import os
import sys
import json
import pprint
import requests
import datetime

""" Download path were all water level time series will be stored """
download_path = 'downloads/'
if not os.path.isdir(download_path):
	os.mkdir(download_path)

args = {}
""" User configuration """
args['username'] = 'username'
args['password'] = 'password'

""" General search options """
args['basin'] = 'Amazon'
#args['continent'] = 'Asia'
#args['country'] = 'de'
#args['min_lon'] = 0
#args['max_lon'] = 10
#args['min_lat'] = 0
#args['max_lat'] = 10

url = 'https://dahiti.dgfi.tum.de/api/v1/'

""" send request as method POST """
args['action'] = 'list-targets'
response = requests.post(url, data=args)
if response.status_code == 200:
	""" convert json string in python list """
	data = json.loads(response.text)
	print ('Dataset(s) found:',len(data))
	i = 0
	for record in data:
		print ('Downloading ... ',record['id'],'->',record['target_name'].encode("utf8"),'('+os.path.abspath(download_path+'/'+str(record['id'])+'.txt')+')')

		""" download water level time series """
		args['action'] = 'download'
		args['dahiti_id'] = record['id']

		response_download = requests.post(url, data=args)
		if response_download.status_code == 200:
			""" convert json string in python list """
			data = json.loads(response_download.text)
			output = open(download_path+'/'+str(record['id'])+'.txt','w')
			output.write('# DAHITI-ID       : '+str(record['id'])+'\n')
			output.write('# Target name     : '+record['target_name'].encode("utf8").decode()+'\n')
			if record['location'] == None:
				output.write('# Location        : '+str(record['location'])+'\n')
			else:
				output.write('# Location        : '+str(record['location'].encode("utf8").decode())+'\n')
			output.write('# Continent       : '+str(record['continent'])+'\n')
			output.write('# Country         : '+str(record['country'])+'\n')
			output.write('# Longitude       : '+str(record['longitude'])+'\n')
			output.write('# Latitude        : '+str(record['latitude'])+'\n')
			output.write('# Points          : '+str(len(data))+'\n')
			output.write('# Software        : '+str(record['water_level']['software'])+'\n')
			output.write('# Download        : '+str(datetime.datetime.now())[0:19]+'\n')
			output.write('# Last-Update     : '+str(record['water_level']['last_update'])+'\n')			
			output.write('# ----------------------------------------\n')
			output.write('# column 1        : date [yyyy-mm-dd]\n')
			output.write('# column 2        : normal heights w.r.t. geoid model [m]\n')
			output.write('# column 3        : error [m] (Kalman Filter, formal errors)\n')
			output.write('# column 4        : dataset (OP:operational, RT:real-time, CB:operational+real-time)\n')
			output.write('# ----------------------------------------\n')
			for entry in data['target']['data']:								
				output.write(str(entry['date'])+' '+str(entry['height'])+' '+str(entry['error'])+' '+str(entry['data_type'])+'\n')
			output.close()
		elif response_download.status_code == 403:
			print ('Status-Code:',response_download.status_code,'- Permission denied!')
		else:
			print ('Status-Code:',response_download.status_code)
		
elif response.status_code == 403:
	print ('Status-Code:',response.status_code,'- Permission denied!')
else:
	print ('Status-Code:',response.status_code)
