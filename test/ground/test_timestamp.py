import unittest
from m4.ground.timestamp import Timestamp


class TestTimestamp(unittest.TestCase):

    def testStringBeginsWith(self):
        ts = Timestamp()
        self.assertTrue(ts.asNowString().startswith(ts.asTodayString(), 0, 8))
        self.assertEqual(str(ts), ts.asNowString())

    # This will fail if now() and today() are called across midnight
    def testStaticMethods(self):
        self.assertTrue(Timestamp.now().startswith(Timestamp.today(), 0, 8))
