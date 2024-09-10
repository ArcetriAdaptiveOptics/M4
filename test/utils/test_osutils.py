'''
Authors
  - C. Selmi: written in 2022
'''
import unittest
import os
from m4.utils.osutils import FileWalker
from pathlib import Path
import tempfile
import shutil


class testFileWalker(unittest.TestCase):

    def setUp(self):
        # creare un albero finto da qualche parte dinamicamente
        self.tmp_path = Path(Path(tempfile.gettempdir()), 'tmp')
        tn_path = Path(self.tmp_path, '20241225_000000')
        os.makedirs(tn_path, exist_ok=True)
        Path(tn_path, '20241225_010000.fits').touch()
        Path(tn_path, '20241225_020000.fits').touch()
        Path(tn_path, '20241226_040000.fits').touch()

        tn_path = Path(self.tmp_path, '20241226_000000')
        os.makedirs(tn_path, exist_ok=True)
        Path(tn_path, '20241226_010000.fits').touch()
        Path(tn_path, '20241226_020000.fits').touch()
        Path(tn_path, '20241226_030000.fits').touch()

        tn_path = Path(self.tmp_path, '20241227_000000')
        os.makedirs(tn_path, exist_ok=True)
        Path(tn_path, '20241227_010000.fits').touch()
        Path(tn_path, '20241227_020000.fits').touch()
        Path(tn_path, '20241227_030000.fits').touch()
        print("Created sample tree into: ", self.tmp_path)

    def testFindTracknum(self):
        # test che il letto da filewalker trova  e fa tutto
        fw = FileWalker(self.tmp_path)
        print("test find tracknum dir 20241227_000000")
        # print("Found tracknum %s" % fw.findTracknum('20241227_000000'))
        self.assertEqual(Path(self.tmp_path, '20241227_000000'),
                         fw.findTracknum('20241227_000000')[0])

        print("test find tracknum file 20241226_020000.fits")
        # print("Found tracknum %s" % fw.findTracknum('20241226_020000'))
        try:
            trovato = fw.findTracknum('20241226_020000')[0]
        except Exception as e:
            print("Exception: ", e)
            raise e
            # trovato = "pippo"

        self.assertEqual(trovato, Path(
            Path(self.tmp_path, '20241226_000000'), '20241226_020000.fits'))

    def testFindTagsBetweenDates(self):
        fw = FileWalker(self.tmp_path)
        print("test find_tags_between_dates 20241226_020000 and  20241227_000000")
        self.assertEqual(
            fw.find_tag_between_dates(
                "20241226_020000", "20241227_000000.fits"),
            [Path(Path(self.tmp_path, '20241225_000000'), '20241226_040000.fits'),
             Path(Path(self.tmp_path, '20241226_000000'), '20241226_020000.fits'),
             Path(Path(self.tmp_path, '20241226_000000'), '20241226_030000.fits'),
             Path(self.tmp_path, '20241227_000000')])

    def tearDown(self):
        try:
            print("Cleaning temporary tree")
            shutil.rmtree(self.tmp_path)
        except Exception as e:
            print("Exception: ", e)
        pass


# if __name__ == "__main__":
#    tf = testFileWalker()
#    tf.setUp()
#    tf.testFindTracknum()
#    tf.testFindTagsBetweenDates()
#    tf.tearDown()
