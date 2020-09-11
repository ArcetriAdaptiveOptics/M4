'''
@author: cselmi
'''

import sys
import time
from opcua import Client
from opcua import ua
from functools import *
from pyzabbix import ZabbixAPI
from pyzabbix import ZabbixMetric, ZabbixSender


def mainZabbix(zserver, port):
    zapi = createZabbixAPI()
    hostname = 'name'
    valore = 7
    valore2 = getZabbixMetrics(zapi, hostname, 'key')

    packet = [
        ZabbixMetric(hostname, 'Descrizione oggetto1', valore),
        # multiple metrics can be sent in same call for effeciency
        ZabbixMetric(hostname, 'Descrizione oggetto2', valore2)
                ]

    result = ZabbixSender(zserver, port, use_config=None).send(packet)
    return result

def createZabbixAPI():
    zapi = ZabbixAPI(url='http://localhost/zabbix/', user='Admin', password='zabbix')
    return zapi

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


#####################
#EndPoint: TcOpcUaServer@CP-4A2AC8[None,None] [opc.tcp://CP-4A2AC8:4840]

def mainOpcUa():
    server = "opc.tcp://192.168.22.100:48050"
    #m4ws = "opc.tcp://10.20.30.6:4840"
    client = Client(url=server)
    client.connect()

    #ok, ma non serve
    root = client.get_root_node()
    print("Root node is: ", root)
    objects = client.get_objects_node()
    print("Objects node is: ", objects)
    print("Children of root are: ", root.get_children())

    #node = root.get_child(["0:Objects", "4:tutti i passaggi dell'albero", "4:fino alla", "4:proprietà che voglio leggere"])
    #var = node.get_value()

    var = client.get_node("ns=7;s=MAIN.i_Temperature_Sensor[16]")
    value = var.get_value()
    type = var.get_data_type_as_variant_type()
    new_value = ua.DataValue(ua.Variant(int(1), type))
    var.set_value(new_value)

#rotazione torre
    rot_angle = client.get_node("ns=7;s=MAIN.f_targetPosition_input[0]")
    type_rot_angle = rot_angle.get_data_type_as_variant_type()
    move_rot = client.get_node("ns=7;s=MAIN.b_MoveCmd[0]")
    tipo_move_rot = move_rot.get_data_type_as_variant_type()

    rot_angle.set_value(1, type_rot_angle)
    move_rot.set_value(True, tipo_move_rot)

    client.disconnect()
    return var
