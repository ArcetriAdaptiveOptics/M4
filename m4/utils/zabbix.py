'''
Created on 4 set 2020

@author: rm
'''


import math
import logging
import sys
logging.basicConfig(stream=sys.stdout, level=logging.INFO)
import time
#from pyzabbix.api import ZabbixAPI
from pyzabbix import ZabbixAPI
from pyzabbix import ZabbixMetric, ZabbixSender
from functools import *
from meteocalc import Temp, dew_point, heat_index, wind_chill, feels_like


def get_frost_point_c(t_air_c, dew_point_c):
    """Compute the frost point in degrees Celsius

    :param t_air_c: current ambient temperature in degrees Celsius
    :type t_air_c: float
    :param dew_point_c: current dew point in degrees Celsius
    :type dew_point_c: float
    :return: the frost point in degrees Celsius
    :rtype: float
    """
    dew_point_k = 273.15 + dew_point_c
    t_air_k = 273.15 + t_air_c
    frost_point_k = dew_point_k - t_air_k + 2671.02 / ((2954.61 / t_air_k) + 2.193665 * math.log(t_air_k) - 13.3448)
    return frost_point_k - 273.15


def get_dew_point_c(t_air_c, rel_humidity):
    """Compute the dew point in degrees Celsius

    :param t_air_c: current ambient temperature in degrees Celsius
    :type t_air_c: float
    :param rel_humidity: relative humidity in %
    :type rel_humidity: float
    :return: the dew point in degrees Celsius
    :rtype: float
    """
    A = 17.27
    B = 237.7
    alpha = ((A * t_air_c) / (B + t_air_c)) + math.log(rel_humidity/100.0)
    return (B * alpha) / (A - alpha)

def getZabbixMetrics(zapi,host, key):
    zversion = zapi.do_request('apiinfo.version')
    print( "Zabbix API version: {}".format(zversion['result']))


    # https://www.zabbix.com/documentation/2.2/manual/api/reference/host/get
    # Get specified host
    print( "----------------------------")
    thehost = zapi.do_request('host.get',
                          {
                              'filter': {'host': host},
                              'selectItems' : 'extend',
                              'output': 'extend'
                          })
    if len(thehost['result'])<1:
        print( "HALTING. There was no host defined in zabbix with id: {}".format(host))
        sys.exit(2)
    hostId = thehost['result'][0]['hostid']
    print( "Found host {} with id {}".format(host,hostId))

    # now look for item within that host
    itemId = None
    for item in thehost['result'][0]['items']:
      # for debugging
      #print( "item[{}] -> {}".format(item['itemid'],item['key_']))
      # if match, then get out int id and type (0=float,1=char,3=unsign,4=text)
        if item['key_'] == key:
            itemId = item['itemid']
            itemType = item['value_type']
    if itemId is None:
        print( "HALTING. There was no item defined on host {} with name: {}".format(host,key))
        sys.exit(2)
    print( "Found item {} on host {} with item id/type {}/{}".format(key,host,itemId,itemType))


    # https://www.zabbix.com/documentation/2.2/manual/api/reference/history/get
    print( "----------------------------")
    history = zapi.do_request('history.get',
                          {
                              'history': itemType,
                              'filter': {'host': host, 'itemid': itemId},
                              'limit': '5',
                              'sortfield': 'clock',
                              'sortorder': 'DESC',
                              'output': 'extend'
                          })
    # show history rows 
    print( "Retrieved {} rows of history".format(len(history['result'])))
    for hist in history['result']:
        # convert epoch to human readable format
        timestr = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(int(hist['clock'])))
        print( "{} @ {}".format(timestr,hist['value']))
    avehist = reduce((lambda x,y: float(x)+float(y)), [d['value'] for d in (history['result'])] ) /len(history['result'])
    print( "Avehist {}.".format(avehist))
    return avehist

def main():
    # Create ZabbixAPI class instance
    zserver = '193.206.155.82'
    port = 10051

    zapi = ZabbixAPI(url='http://localhost/zabbix/', user='Admin', password='zabbix')

    # Get all monitored hosts
    hostname2 = 'moxa_host2'
    hostname3 = 'moxa_host3'
    h0 = getZabbixMetrics(zapi, hostname3,'IHA-ERIS-HUM')
    h1 = getZabbixMetrics(zapi, hostname3,'IHA-HYDAC-HUM')
    t0 = getZabbixMetrics(zapi, hostname2,'IHA-DSMSIM-Board-RTD')
    t1 = getZabbixMetrics(zapi, hostname2,'IHA-HYDAC-RTD')
    print("ERIS T {} and H {}, HYDAC T {} and H {}".format(t0, h0, t1, h1))


    DewPointERIS  = dew_point(Temp(t0,'c'), humidity=h0)
    DewPointHYDAC = dew_point(Temp(t1,'c'), humidity=h1)
    print("Dew Point ERIS {} and HYDAC {}".format(DewPointERIS,DewPointHYDAC))
    # Send metrics to zabbix trapper
    packet = [
        ZabbixMetric(hostname3, 'DewPointERIS',DewPointERIS)
        # multiple metrics can be sent in same call for effeciency
        ,ZabbixMetric(hostname3, 'DewPointHYDAC', DewPointHYDAC)
    ]

    result = ZabbixSender(zserver,port,use_config=None).send(packet)
    print(result)



main()