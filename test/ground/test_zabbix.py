'''
Authors
  - C. Selmi: written in 2021
'''
import unittest
from unittest.mock import patch


class Test(unittest.TestCase):

    @patch('opcua.Client', autospec=True)
    @patch('pyzabbix.ZabbixAPI', autospec=True)
    @patch('pyzabbix.ZabbixSender', autospec=True)
#    @patch('pyzabbix.ZabbixMetric', autospec=True)
    def setUp(self, mock_client, mock_api, mock_sender):
        from m4.ground import zabbix
        self.zabbix = zabbix


    def tearDown(self):
        self.zabbix

    @unittest.skip(" module pyzabbix from /opt/hostedtoolcache/Python/3.7.10/x64/lib/python3.7/site-packages/pyzabbix/__init__.pydoes not have the attribute ZabbixSender ")
    def testZabbix(self):
        self.zabbix.mainZabbix()
