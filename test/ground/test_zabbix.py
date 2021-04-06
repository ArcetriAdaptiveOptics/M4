'''
Authors
  - C. Selmi: written in 2021
'''
import unittest
import mock


class Test(unittest.TestCase):

    @mock.patch('opcua.Client', autospec=True)
    @mock.patch('pyzabbix.ZabbixAPI', autospec=True)
    @mock.patch('pyzabbix.ZabbixSender', autospec=True)
#    @mock.patch('pyzabbix.ZabbixMetric', autospec=True)
    def setUp(self, mock_client, mock_api, mock_sender):
        from m4.ground import zabbix
        self.zabbix = zabbix


    def tearDown(self):
        self.zabbix

    @unittest.skip(" module pyzabbix from /opt/hostedtoolcache/Python/3.7.10/x64/lib/python3.7/site-packages/pyzabbix/__init__.pydoes not have the attribute ZabbixSender ")
    def testZabbix(self):
        self.zabbix.mainZabbix()
