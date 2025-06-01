'''
Authors
  - C. Selmi: written in 2021
'''
import unittest
from m4.simulator.fake_temperature_sensors import FakeTemperatureSensors

class TestFakeTemperatureSensor(unittest.TestCase):

    def testGetTemperature(self):
        fake = FakeTemperatureSensors()
        temp = fake.getTemperature()
        message = "Get position doesn't return a vector of shape 24."
        self.assertEqual(temp.shape, (24,), message)
