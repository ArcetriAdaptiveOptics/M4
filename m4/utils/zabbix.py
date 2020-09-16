'''
@author: cselmi
'''

import sys
import time
import numpy as np
from functools import *
from pyzabbix import ZabbixAPI
from pyzabbix import ZabbixMetric, ZabbixSender


def mainZabbix():
    zapi = createZabbixAPI()
    hostname = 'M4OTT'
    zserver = '192.168.22.22'
    port = 10051

    valore2 = getZabbixMetrics(zapi, hostname, 'key')

    packet = [
        ZabbixMetric(hostname, 'Nome oggetto1', valore),
        # multiple metrics can be sent in same call for effeciency
        ZabbixMetric(hostname, 'Nome oggetto2', valore2)
                ]

    result = ZabbixSender(zserver, port, use_config=None).send(packet)
    return result

def createZabbixAPI():
    zapi = ZabbixAPI(url='http://192.168.22.22/zabbix/', user='Admin', password='zabbix')
    return zapi

def read_temperature(zapi, hostname):
    temp_list = []
    for i in range(24):
        temp = getZabbixMetrics(zapi, hostname, 'PT%02d' %i)
        temp_list.append(temp)
    return np.array(temp_list)

def getZabbixMetrics(zapi, host, key):
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
                              'filter': {'itemid': itemId},
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
