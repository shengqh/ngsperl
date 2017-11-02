"""@package Oasis API
Documentation for Oasis smallRNA API.

It launches one or several zip files to oasis.dzne.de to upload 
the files and assess the quality of the uploaded sequencing data and summarize the
 quantity of reads for each small RNA species in one count file per samples.
____ ____ ____ _ ____    ____ ___  _ 
|  | |__| [__  | [__     |__| |__] | 
|__| |  | ___] | ___]    |  | |    | 
                                     


INTRODUCTION
------------

This API is developed to allow multiple batches submission to the Oasis 
web server. It works using Python.

Copyright (C) 2015  Bonnlab

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import os
import sys
import glob
try:
	import zipfile
except ImportError:
	raise ImportError('Please install zipfile module')
try:
	import requests
except ImportError:
	raise ImportError('Please install requests module')
import time
try:
	import requests_toolbelt
except ImportError:
	raise ImportError('Please install requests_toolbelt module')
from requests_toolbelt import MultipartEncoder, MultipartEncoderMonitor

parser = argparse.ArgumentParser(description='Send job in VSAC Web API.', \
epilog="Welcome to OASIS. \
Launching this part of the analysis suite will allow you to assess the quality \
of your sequencing experiment and summarize the quantity of reads for each \
small RNA species in one count file per sample. These count files can then be \
used to perform differential expression analysis or biomarker discovery. ")
parser.add_argument("-g", "--genome", help="Organism your samples originate \
from. Values can be: bostau8 for Bos taurus, canfam3 for Canis familiaris, \
ce10 for Caenorhabditis elegans, danrer10 for Danio rerio, dm6 for Drosophila \
melanogaster, equcab2 for Equus caballus, galgal4 for Gallus gallus, gorgor3 for Gorilla gorilla, hg38 for\
 Homo sapiens, mm10 for Mus musculus, pantro4 for Pan troglodytes, rn6 for \
Rattus norvegicus and ss3 for Sus scrofa.", required=True, \
dest="genome", type=str, default="hs", action='store', choices=['bostau8', 'canfam3', 'ce10', \
'danrer10', 'dm6', 'equcab2','galgal4', 'gorgor3', 'hg38', 'mm10', 'pantro4', 'rn6', 'ss3'])
parser.add_argument("-a", "--adapter", help="Adapter that was used to prepare \
the libraries.", required=True)
parser.add_argument("-e", "--email", help="Email to send output of analysis.", required=True)
parser.add_argument("-s", "--samples", nargs='+', help="Zipped directory with the \
samples created with our tool FastQcCompressor.jar", required=True)
parser.add_argument("-en", "--experimentName", help="A string defining the name \
of the experiment", required=False, type=str, default='')
parser.add_argument("-m", "--mismatches", help="Please provide the percentage of the read length to use as number of \
mismatches in the alignment phase", required=False, type=str, action='store', default='5', choices=['0', '5', '10', '20'])
parser.add_argument("-min", "--minimum", help="Please provide the minimum length \
of the read to use in the analysis", required=False, type=str, default='15', action='store', choices=['15', '16', '17', '18'])
parser.add_argument("-max", "--maximum", help="Please provide the maximum length \
of the read to use in the analysis", required=False, type=str, default='32', action='store', choices=map(str,list(range(32,50))))

requests.packages.urllib3.disable_warnings()
# API_URL='/apismallrnapipelinewww.php'
API_URL='/apismallrnapipelinewww.php'

try:
	args = parser.parse_args()
except argparse.ArgumentError:
	print("Unexpected error:", sys.exc_info()[0])
	sys.exit(1)
version='1.1.0'
try:
	r=requests.get('https://oasis.dzne.de/apiversion.php', verify=False)
except:
	print("Oasis is unavailable, please contact to oasis@dzne.de. Unexpected error:", sys.exc_info()[0])
	sys.exit(1)

if r.text!=version:
	print("Please download the latest OASIS API version from https://oasis.dzne.de")
	sys.exit(1)
extensions=['.dat', '.stats', '.count', '.fastq', '.fastq.gz','.zip']
if len(args.samples)==1:
	f=''.join(args.samples)
	m = MultipartEncoder(fields={'reference': args.genome, 'adapter': args.adapter, \
	'experimentName':args.experimentName, 'mismatch':args.mismatches, 'minread':args.minimum, \
	'maxread':args.maximum, 'email':args.email, 'uploaded_file': ('uploaded_file', open(f, 'rb'), 'application/zip')})
	print ("It's uploading..."+ time.strftime("%D %H:%M"))
	r=requests.post('https://oasis.dzne.de'+API_URL,data=m, headers={'Content-Type': m.content_type}, verify=False)
	try:
		if r.status_code==requests.codes.ok:
			print ("Uploaded at " + time.strftime("%D %H:%M"))
			print (r.text)
	except:
		print("Unexpected error ESTE:", sys.exc_info()[0])
		sys.exit(1)
	print ("Process zip file: " + ''.join(args.samples))
	sys.exit(0)
elif len(args.samples)>1:
	"""Batch process of several zip files"""
	for f in args.samples:
		m = MultipartEncoder(fields={'reference': args.genome, 'adapter': args.adapter,\
		'experimentName':args.experimentName, 'mismatch':args.mismatches, 'minread':args.minimum, \
		'maxread':args.maximum, 'email':args.email, 'uploaded_file': ('uploaded_file', open(f, 'rb'), 'application/zip')})
		files={'uploaded_file':open(f, 'rb')}
		print ("Proceed to upload " + f)
		try:
			print ("It's uploading..."+ time.strftime("%D %H:%M"))
			r=requests.post('https://oasis.dzne.de'+API_URL,data=m, headers={'Content-Type': m.content_type, 'User-agent':'smallRNA.py/0.4'})
			if r.status_code==requests.codes.ok:
				print ("Uploaded at " + time.strftime("%D %H:%M"))
				print (r.text)
		except:
			print("Unexpected error ESTEOTRO:", sys.exc_info()[0])
			sys.exit(1)
sys.exit(0)

