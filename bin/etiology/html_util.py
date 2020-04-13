#! /usr/bin/python3
# -*- coding: utf-8 -*-

import base64
import os


class HtmlUtil:
	
	@staticmethod
	def getPNGBinary(name):
	    encoded_string = ''
	    with open(os.path.dirname(os.path.realpath(__file__)) + '/html_template/images/' + name, 'rb') as image_file:
	        encoded_string = base64.b64encode(image_file.read())
	    # print("data:image/png;base64," + str(encoded_string)[2:-1])
	    return "data:image/png;base64," + str(encoded_string)[2:-1]

	@staticmethod    
	def getFileContent(name):
	    strContent = ''
	    with open(os.path.dirname(os.path.realpath(__file__)) + '/html_template/css/'+name, 'r') as file:
	        for line in file.readlines():
	            strContent += line
	    return strContent