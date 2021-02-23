'''
Authors
  - C. Selmi: written in 2021
'''
import unittest
import mock


class Test(unittest.TestCase):

    @mock.patch('opcua.Client', autospec=True)
    def setUp(self, mock_client):
        from m4.ground import zabbix
        self.zabbix = zabbix


    def tearDown(self):
        pass

#    @unittest.skip('Non funziona')
    @mock.patch('pyzabbix.ZabbixAPI', autospec=True)
    @mock.patch('pyzabbix.ZabbixSender', autospec=True)
    @mock.patch('pyzabbix.ZabbixMetric', autospec=True)
    def testZabbix(self, mock_api, mock_sender, mock_metric):
        self.zabbix.mainZabbix()
