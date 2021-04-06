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
    @mock.patch('pyzabbix.ZabbixMetric', autospec=True)
    def setUp(self, mock_client, mock_api, mock_sender, mock_metric):
        from m4.ground import zabbix
        self.zabbix = zabbix


    def tearDown(self):
        self.zabbix

    def testZabbix(self):
        self.zabbix.mainZabbix()
